import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import shapiro, ttest_ind, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

print("=== METRIKA-MOCA STATISZTIKAI ELEMZÉS ===")

# 1. RANGES betöltése
print("\n1. Ranges betöltése...")
# 1. RANGES betöltése
print("\n1. Ranges betöltése...")
try:
    # Egyszerűen pandas-szal olvassuk be - vesszővel elválasztott CSV
    ranges_df = pd.read_csv('ranges.csv', header=None, names=['Metric', 'OptimalMin', 'OptimalMax', 'ValidMin', 'ValidMax'])
    
    print(f"Ranges betöltve: {len(ranges_df)} metrika")
    
    if len(ranges_df) > 0:
        print(f"Első pár metrika:")
        for i in range(min(5, len(ranges_df))):
            row = ranges_df.iloc[i]
            print(f"  {row['Metric']}: opt[{row['OptimalMin']}-{row['OptimalMax']}], valid[{row['ValidMin']}-{row['ValidMax']}]")
    else:
        print("HIBA: Egyetlen metrika sem lett betöltve!")
        exit(1)
    
except Exception as e:
    print(f"HIBA ranges betöltésnél: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
    
  

# 2. Fő adatok betöltése
print("\n2. Főadatok betöltése...")
try:
    # Metrikák betöltése
    data_df = pd.read_csv('specific_indices_matrix.csv')
    print(f"Metrika adatok betöltve: {len(data_df)} páciens")
    
    # MOCA adatok betöltése
    moca_df = pd.read_csv('data/df_unnormalized.csv')
    print(f"MOCA fájl betöltve: {len(moca_df)} sor")
    
    # Identifier mentén unique-vá tétel
    moca_unique = moca_df.drop_duplicates(subset=['Identifier'])
    print(f"Unique páciensek MOCA adatok: {len(moca_unique)}")
    
    # DeltaMOCA1 oszlop kiválasztása
    if 'DeltaMOCA1' not in moca_unique.columns:
        print("HIBA: DeltaMOCA1 oszlop nem található!")
        print(f"Elérhető oszlopok: {list(moca_unique.columns)}")
        exit(1)
    
    moca_selected = moca_unique[['Identifier', 'DeltaMOCA1']].copy()
    print(f"MOCA adatok kiválasztva: {len(moca_selected)} páciens")
    
    # Merge a metrika adatokkal
    data_df = data_df.merge(moca_selected, left_on='PatientID', right_on='Identifier', how='inner')
    print(f"Merge után: {len(data_df)} páciens")
    
    moca_col = 'DeltaMOCA1'
    print(f"MOCA oszlop: {moca_col}")
    
    # MOCA adatok ellenőrzése
    valid_moca = data_df[moca_col].notna()
    print(f"Valid MOCA értékek: {sum(valid_moca)}/{len(data_df)}")
    
except Exception as e:
    print(f"HIBA adatok betöltésnél: {e}")
    exit(1)

# 3. Eredmény tárolók
results_metric_to_moca = []  # Metrika jó/rossz → MOCA különbség
results_moca_to_metric = []  # MOCA jó/rossz → Metrika különbség

print(f"\n3. Elemzés kezdése {len(ranges_df)} metrikára...")

# 4. FŐCIKLUS - Minden metrikán végigmegy
for idx, range_row in ranges_df.iterrows():
    metric_name = range_row['Metric']
    opt_min = range_row['OptimalMin']
    opt_max = range_row['OptimalMax']
    valid_min = range_row['ValidMin']
    valid_max = range_row['ValidMax']
    
    print(f"\n--- {idx+1}/{len(ranges_df)}: {metric_name} ---")
    
    # Metrika oszlop keresése az adatokban
    if metric_name not in data_df.columns:
        print(f"  SKIP: Metrika nem található az adatokban")
        continue
    
    # Valid adatok szűrése
    mask_valid = (
        data_df[metric_name].notna() & 
        data_df[moca_col].notna() &
        (data_df[metric_name] >= valid_min) & 
        (data_df[metric_name] <= valid_max)
    )
    
    valid_data = data_df[mask_valid].copy()
    
    if len(valid_data) < 6:  # Minimum 6 adat kell (3-3 csoportonként)
        print(f"  SKIP: Túl kevés valid adat ({len(valid_data)})")
        continue
    
    print(f"  Valid adatok: {len(valid_data)}")
    
    # =================================================================
    # ELEMZÉS 1: METRIKA JÓ/ROSSZ → MOCA KÜLÖNBSÉG
    # =================================================================
    
    # Metrika alapján csoportok: optimális tartományban vs kívül
    mask_good_metric = (
        (valid_data[metric_name] >= opt_min) & 
        (valid_data[metric_name] <= opt_max)
    )
    
    good_metric_group = valid_data[mask_good_metric][moca_col].values
    bad_metric_group = valid_data[~mask_good_metric][moca_col].values
    
    if len(good_metric_group) >= 3 and len(bad_metric_group) >= 3:
        
        # Normalitás tesztek
        _, good_p_norm = shapiro(good_metric_group) if len(good_metric_group) <= 5000 else (None, 0.001)
        _, bad_p_norm = shapiro(bad_metric_group) if len(bad_metric_group) <= 5000 else (None, 0.001)
        
        good_normal = good_p_norm > 0.05 if good_p_norm is not None else False
        bad_normal = bad_p_norm > 0.05 if bad_p_norm is not None else False
        
        # Teszt választása és végrehajtása
        if good_normal and bad_normal:
            # Mindkettő normális → t-test
            stat, p_value = ttest_ind(good_metric_group, bad_metric_group)
            test_type = "t-test"
            
            # Cohen's d számítása
            pooled_std = np.sqrt(((len(good_metric_group)-1)*np.var(good_metric_group, ddof=1) + 
                                 (len(bad_metric_group)-1)*np.var(bad_metric_group, ddof=1)) / 
                                (len(good_metric_group) + len(bad_metric_group) - 2))
            cohens_d = (np.mean(good_metric_group) - np.mean(bad_metric_group)) / pooled_std if pooled_std > 0 else 0
            effect_size = abs(cohens_d)
            
        else:
            # Legalább egy nem normális → Mann-Whitney U
            stat, p_value = mannwhitneyu(good_metric_group, bad_metric_group, alternative='two-sided')
            test_type = "Mann-Whitney"
            
            # r effect size számítása
            z_score = stats.norm.ppf(p_value/2)  # közelítés
            effect_size = abs(z_score) / np.sqrt(len(good_metric_group) + len(bad_metric_group))
        
        # Szignifikancia
        significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        
        # Eredmény tárolása
        results_metric_to_moca.append({
            'Metric': metric_name,
            'Analysis': 'Metric_to_MOCA',
            'Test_Type': test_type,
            'p_value': p_value,
            'Effect_Size': effect_size,
            'Significance': significance,
            'Good_Metric_N': len(good_metric_group),
            'Good_Metric_MOCA_Mean': np.mean(good_metric_group),
            'Good_Metric_MOCA_Std': np.std(good_metric_group, ddof=1),
            'Bad_Metric_N': len(bad_metric_group),
            'Bad_Metric_MOCA_Mean': np.mean(bad_metric_group),
            'Bad_Metric_MOCA_Std': np.std(bad_metric_group, ddof=1),
            'Good_Normal': good_normal,
            'Bad_Normal': bad_normal
        })
        
        print(f"    Metrika→MOCA: {test_type}, p={p_value:.4f}{significance}, effect={effect_size:.3f}")
        print(f"      Jó metrika MOCA: μ={np.mean(good_metric_group):.3f} (n={len(good_metric_group)})")
        print(f"      Rossz metrika MOCA: μ={np.mean(bad_metric_group):.3f} (n={len(bad_metric_group)})")
    
    else:
        print(f"    Metrika→MOCA: SKIP (jó n={len(good_metric_group)}, rossz n={len(bad_metric_group)})")
    
    # =================================================================
    # ELEMZÉS 2: MOCA JÓ/ROSSZ → METRIKA KÜLÖNBSÉG
    # =================================================================
    
    # MOCA alapján csoportok
    mask_good_moca = valid_data[moca_col] >= -2
    
    good_moca_group = valid_data[mask_good_moca][metric_name].values
    bad_moca_group = valid_data[~mask_good_moca][metric_name].values
    
    if len(good_moca_group) >= 3 and len(bad_moca_group) >= 3:
        
        # Normalitás tesztek
        _, good_p_norm = shapiro(good_moca_group) if len(good_moca_group) <= 5000 else (None, 0.001)
        _, bad_p_norm = shapiro(bad_moca_group) if len(bad_moca_group) <= 5000 else (None, 0.001)
        
        good_normal = good_p_norm > 0.05 if good_p_norm is not None else False
        bad_normal = bad_p_norm > 0.05 if bad_p_norm is not None else False
        
        # Teszt választása és végrehajtása
        if good_normal and bad_normal:
            # Mindkettő normális → t-test
            stat, p_value = ttest_ind(good_moca_group, bad_moca_group)
            test_type = "t-test"
            
            # Cohen's d számítása
            pooled_std = np.sqrt(((len(good_moca_group)-1)*np.var(good_moca_group, ddof=1) + 
                                 (len(bad_moca_group)-1)*np.var(bad_moca_group, ddof=1)) / 
                                (len(good_moca_group) + len(bad_moca_group) - 2))
            cohens_d = (np.mean(good_moca_group) - np.mean(bad_moca_group)) / pooled_std if pooled_std > 0 else 0
            effect_size = abs(cohens_d)
            
        else:
            # Legalább egy nem normális → Mann-Whitney U
            stat, p_value = mannwhitneyu(good_moca_group, bad_moca_group, alternative='two-sided')
            test_type = "Mann-Whitney"
            
            # r effect size számítása
            z_score = stats.norm.ppf(p_value/2) if p_value > 0 else 3  # közelítés
            effect_size = abs(z_score) / np.sqrt(len(good_moca_group) + len(bad_moca_group))
        
        # Szignifikancia
        significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        
        # Eredmény tárolása
        results_moca_to_metric.append({
            'Metric': metric_name,
            'Analysis': 'MOCA_to_Metric',
            'Test_Type': test_type,
            'p_value': p_value,
            'Effect_Size': effect_size,
            'Significance': significance,
            'Good_MOCA_N': len(good_moca_group),
            'Good_MOCA_Metric_Mean': np.mean(good_moca_group),
            'Good_MOCA_Metric_Std': np.std(good_moca_group, ddof=1),
            'Bad_MOCA_N': len(bad_moca_group),
            'Bad_MOCA_Metric_Mean': np.mean(bad_moca_group),
            'Bad_MOCA_Metric_Std': np.std(bad_moca_group, ddof=1),
            'Good_Normal': good_normal,
            'Bad_Normal': bad_normal
        })
        
        print(f"    MOCA→Metrika: {test_type}, p={p_value:.4f}{significance}, effect={effect_size:.3f}")
        print(f"      Jó MOCA metrika: μ={np.mean(good_moca_group):.3f} (n={len(good_moca_group)})")
        print(f"      Rossz MOCA metrika: μ={np.mean(bad_moca_group):.3f} (n={len(bad_moca_group)})")
    
    else:
        print(f"    MOCA→Metrika: SKIP (jó MOCA n={len(good_moca_group)}, rossz MOCA n={len(bad_moca_group)})")

# =================================================================
# EREDMÉNYEK MENTÉSE ÉS ÖSSZEFOGLALÁS
# =================================================================

print(f"\n=== EREDMÉNYEK ÖSSZEFOGLALÁSA ===")

# DataFrame-ek létrehozása
df_metric_to_moca = pd.DataFrame(results_metric_to_moca)
df_moca_to_metric = pd.DataFrame(results_moca_to_metric)

print(f"Metrika→MOCA elemzések: {len(df_metric_to_moca)}")
print(f"MOCA→Metrika elemzések: {len(df_moca_to_metric)}")

# Szignifikáns eredmények
if len(df_metric_to_moca) > 0:
    sig_metric_to_moca = df_metric_to_moca[df_metric_to_moca['p_value'] < 0.05]
    print(f"Szignifikáns Metrika→MOCA: {len(sig_metric_to_moca)}")
    
    if len(sig_metric_to_moca) > 0:
        print("\nTop szignifikáns Metrika→MOCA eredmények:")
        top_results = sig_metric_to_moca.nsmallest(5, 'p_value')
        for _, row in top_results.iterrows():
            print(f"  {row['Metric']}: p={row['p_value']:.4f}{row['Significance']}, effect={row['Effect_Size']:.3f}")

if len(df_moca_to_metric) > 0:
    sig_moca_to_metric = df_moca_to_metric[df_moca_to_metric['p_value'] < 0.05]
    print(f"Szignifikáns MOCA→Metrika: {len(sig_moca_to_metric)}")
    
    if len(sig_moca_to_metric) > 0:
        print("\nTop szignifikáns MOCA→Metrika eredmények:")
        top_results = sig_moca_to_metric.nsmallest(5, 'p_value')
        for _, row in top_results.iterrows():
            print(f"  {row['Metric']}: p={row['p_value']:.4f}{row['Significance']}, effect={row['Effect_Size']:.3f}")

# Fájlba mentés
import os
current_dir = os.getcwd()
print(f"\n=== FÁJLBA MENTÉS ===")
print(f"Jelenlegi könyvtár: {current_dir}")

if len(df_metric_to_moca) > 0:
    filename = 'metric_to_moca_analysis.csv'
    df_metric_to_moca.to_csv(filename, index=False)
    filepath = os.path.join(current_dir, filename)
    print(f"Metrika→MOCA: {filepath} ({len(df_metric_to_moca)} eredmény)")
    print(f"  Fájl létezik: {os.path.exists(filepath)}")
else:
    print("Metrika→MOCA: NINCS ADAT - fájl nem létrehozva")

if len(df_moca_to_metric) > 0:
    filename = 'moca_to_metric_analysis.csv'
    df_moca_to_metric.to_csv(filename, index=False)
    filepath = os.path.join(current_dir, filename)
    print(f"MOCA→Metrika: {filepath} ({len(df_moca_to_metric)} eredmény)")
    print(f"  Fájl létezik: {os.path.exists(filepath)}")
else:
    print("MOCA→Metrika: NINCS ADAT - fájl nem létrehozva")

# Kombinált eredmények
if len(df_metric_to_moca) > 0 and len(df_moca_to_metric) > 0:
    filename = 'combined_metric_moca_analysis.csv'
    combined_results = pd.concat([df_metric_to_moca, df_moca_to_metric], ignore_index=True)
    combined_results.to_csv(filename, index=False)
    filepath = os.path.join(current_dir, filename)
    print(f"Kombinált: {filepath} ({len(combined_results)} eredmény)")
    print(f"  Fájl létezik: {os.path.exists(filepath)}")

# Szignifikáns eredmények külön
if len(df_metric_to_moca) > 0 or len(df_moca_to_metric) > 0:
    sig_results = []
    if len(df_metric_to_moca) > 0:
        sig_results.append(df_metric_to_moca[df_metric_to_moca['p_value'] < 0.05])
    if len(df_moca_to_metric) > 0:
        sig_results.append(df_moca_to_metric[df_moca_to_metric['p_value'] < 0.05])
    
    if len(sig_results) > 0:
        all_significant = pd.concat([df for df in sig_results if len(df) > 0], ignore_index=True)
        
        if len(all_significant) > 0:
            filename = 'significant_metric_moca_results.csv'
            all_significant.to_csv(filename, index=False)
            filepath = os.path.join(current_dir, filename)
            print(f"Szignifikáns: {filepath} ({len(all_significant)} eredmény)")
            print(f"  Fájl létezik: {os.path.exists(filepath)}")
        else:
            print("Szignifikáns: NINCS SZIGNIFIKÁNS EREDMÉNY - fájl nem létrehozva")
    else:
        print("Szignifikáns: NINCS ADAT - fájl nem létrehozva")

print(f"\n=== ELEMZÉS BEFEJEZVE ===")
print("Kimeneti fájlok:")
print("  - metric_to_moca_analysis.csv")
print("  - moca_to_metric_analysis.csv") 
print("  - combined_metric_moca_analysis.csv")
print("  - significant_metric_moca_results.csv")
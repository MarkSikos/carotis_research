import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import (
    shapiro, ttest_ind, mannwhitneyu, levene, bartlett, 
    chi2_contingency, fisher_exact, pearsonr, spearmanr
)
import warnings
warnings.filterwarnings('ignore')

print("=== BAND ÉS SCALE-SPECIFIKUS BIOMARKER-MOCA STATISZTIKAI ELEMZÉS ===")

# 1. ADATOK BETÖLTÉSE
print("\n1. Adatok betöltése...")
try:
    # Comprehensive biomarker adatok
    biomarker_df = pd.read_csv('comprehensive_biomarkers_and_stats.csv')
    print(f"Biomarker adatok betöltve: {len(biomarker_df)} sor")
    print(f"Oszlopok: {len(biomarker_df.columns)}")
    
    # MOCA adatok betöltése
    moca_df = pd.read_csv('data/df_unnormalized.csv')
    print(f"MOCA fájl betöltve: {len(moca_df)} sor")
    
    # Unique páciensek
    moca_unique = moca_df.drop_duplicates(subset=['Identifier'])
    print(f"Unique páciensek MOCA adatok: {len(moca_unique)}")
    
    # DeltaMOCA1 oszlop ellenőrzése
    if 'DeltaMOCA1' not in moca_unique.columns:
        print("HIBA: DeltaMOCA1 oszlop nem található!")
        print(f"Elérhető oszlopok: {list(moca_unique.columns)}")
        exit(1)
    
    moca_selected = moca_unique[['Identifier', 'DeltaMOCA1']].copy()
    
    # Unique scale-ek és band-ek azonosítása
    unique_scales = biomarker_df['Scale'].unique()
    unique_bands = biomarker_df['Band'].unique()
    print(f"Elérhető scale-ek: {list(unique_scales)}")
    print(f"Elérhető band-ek: {list(unique_bands)}")
    
    # Kombinációk száma
    total_combinations = len(unique_scales) * len(unique_bands)
    print(f"Összes Band×Scale kombináció: {total_combinations}")
    
except Exception as e:
    print(f"HIBA adatok betöltésnél: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# 2. METRIKA TÍPUSFELISMERÉS ÉS TESZT VÁLASZTÁS
def detect_metric_type(column_name):
    """Metrika típusának felismerése az oszlopnév alapján"""
    
    col_lower = column_name.lower()
    
    # 1. Középérték típusú metrikák
    if any(x in col_lower for x in ['_mean', '_median', '_mode_angle', '_q1', '_q3', '_min', '_max']):
        return 'location'
    
    # 2. Szórás/variancia típusú metrikák  
    if any(x in col_lower for x in ['_std', '_iqr', '_range']):
        return 'variability'
    
    # 3. Eloszlás alakja
    if any(x in col_lower for x in ['_kurtosis', '_skewness']):
        return 'distribution_shape'
    
    # 4. Arány/hányados típusú
    if any(x in col_lower for x in ['_ratio', '_asymmetry']):
        return 'proportion'
    
    # 5. Entrópia és komplexitás
    if any(x in col_lower for x in ['_entropy', '_n_peaks']):
        return 'complexity'
    
    # 6. Korreláció típusú
    if any(x in col_lower for x in ['correlation', '_lag']):
        return 'correlation'
    
    # 7. Spektrális teljesítmény
    if any(x in col_lower for x in ['_power', 'band_power']):
        return 'power'
    
    # 8. Biomarker specifikus
    if any(x in col_lower for x in ['regularity', 'dysfunction', 'capacity', 'activity', 'sensitivity', 'integrity', 'stiffness']):
        return 'biomarker'
    
    # 9. Minőségi/technikai
    if any(x in col_lower for x in ['actualfs', 'dataquality', 'signallength']):
        return 'technical'
    
    # 10. Default
    return 'general'

def choose_statistical_test(metric_type, good_group, bad_group):
    """Statisztikai teszt választása metrika típus alapján"""
    
    # Alapvető ellenőrzések
    if len(good_group) < 3 or len(bad_group) < 3:
        return None, None, None, None, None
    
    # NaN értékek eltávolítása
    good_clean = good_group.dropna()
    bad_clean = bad_group.dropna()
    
    if len(good_clean) < 3 or len(bad_clean) < 3:
        return None, None, None, None, None
    
    try:
        if metric_type == 'location':
            # Középérték típusú metrikák - location tesztek
            return perform_location_test(good_clean, bad_clean)
            
        elif metric_type == 'variability':
            # Szórás típusú metrikák - variability tesztek
            return perform_variability_test(good_clean, bad_clean)
            
        elif metric_type == 'distribution_shape':
            # Eloszlás alakja - shape tesztek
            return perform_location_test(good_clean, bad_clean)  # Kurtosis/Skewness -> location test
            
        elif metric_type == 'proportion':
            # Arány típusú - proportion tesztek
            return perform_proportion_test(good_clean, bad_clean)
            
        elif metric_type == 'complexity':
            # Komplexitás (entropy, peaks) - általában location
            return perform_location_test(good_clean, bad_clean)
            
        elif metric_type == 'correlation':
            # Korreláció - correlation tesztek
            return perform_correlation_test(good_clean, bad_clean)
            
        elif metric_type == 'power':
            # Spektrális teljesítmény - power tesztek
            return perform_power_test(good_clean, bad_clean)
            
        elif metric_type == 'biomarker':
            # Biomarker - általában location
            return perform_location_test(good_clean, bad_clean)
            
        elif metric_type == 'technical':
            # Technikai metrikák - location
            return perform_location_test(good_clean, bad_clean)
            
        else:  # 'general'
            # Általános - location test
            return perform_location_test(good_clean, bad_clean)
            
    except Exception as e:
        print(f"      Teszt hiba: {e}")
        return None, None, None, None, None

def perform_location_test(good_data, bad_data):
    """Központi tendencia tesztek (t-test vagy Mann-Whitney U)"""
    
    # Normalitás tesztek
    good_normal = check_normality(good_data)
    bad_normal = check_normality(bad_data)
    
    if good_normal and bad_normal:
        # Mindkettő normális -> t-test
        stat, p_value = ttest_ind(good_data, bad_data)
        test_type = "t-test"
        
        # Cohen's d effect size
        pooled_std = np.sqrt(((len(good_data)-1)*np.var(good_data, ddof=1) + 
                             (len(bad_data)-1)*np.var(bad_data, ddof=1)) / 
                            (len(good_data) + len(bad_data) - 2))
        effect_size = abs(np.mean(good_data) - np.mean(bad_data)) / pooled_std if pooled_std > 0 else 0
        
    else:
        # Legalább egy nem normális -> Mann-Whitney U
        stat, p_value = mannwhitneyu(good_data, bad_data, alternative='two-sided')
        test_type = "Mann-Whitney U"
        
        # r effect size (rank correlation)
        z_score = abs(stats.norm.ppf(p_value/2)) if p_value > 0 else 3
        effect_size = z_score / np.sqrt(len(good_data) + len(bad_data))
    
    return test_type, p_value, effect_size, good_normal and bad_normal, f"Good normal: {good_normal}, Bad normal: {bad_normal}"

def perform_variability_test(good_data, bad_data):
    """Variancia/szórás tesztek (Levene test)"""
    
    # Levene test (robusztusabb mint F-test)
    stat, p_value = levene(good_data, bad_data)
    test_type = "Levene test"
    
    # Effect size: ratio of variances
    good_var = np.var(good_data, ddof=1)
    bad_var = np.var(bad_data, ddof=1)
    effect_size = max(good_var, bad_var) / min(good_var, bad_var) if min(good_var, bad_var) > 0 else 1
    
    return test_type, p_value, effect_size, True, f"Variance ratio: {effect_size:.3f}"

def perform_proportion_test(good_data, bad_data):
    """Arány/proportion tesztek"""
    
    # Ha értékek 0-1 között vannak, proportion test
    if all(0 <= x <= 1 for x in good_data) and all(0 <= x <= 1 for x in bad_data):
        # Proportion test - használjuk mint location test
        return perform_location_test(good_data, bad_data)
    else:
        # Általános location test
        return perform_location_test(good_data, bad_data)

def perform_correlation_test(good_data, bad_data):
    """Korreláció tesztek"""
    
    # Korrelációs értékekre Fisher Z-transform alapú teszt
    # De gyakorlatilag location test
    return perform_location_test(good_data, bad_data)

def perform_power_test(good_data, bad_data):
    """Spektrális teljesítmény tesztek"""
    
    # Power értékek általában log-normal eloszlásúak
    # Log transform után location test
    try:
        # Log transform (csak pozitív értékekre)
        good_log = np.log(good_data[good_data > 0])
        bad_log = np.log(bad_data[bad_data > 0])
        
        if len(good_log) >= 3 and len(bad_log) >= 3:
            return perform_location_test(good_log, bad_log)
        else:
            return perform_location_test(good_data, bad_data)
    except:
        return perform_location_test(good_data, bad_data)

def check_normality(data, alpha=0.05):
    """Normalitás ellenőrzése"""
    if len(data) <= 3:
        return False
    
    try:
        _, p_value = shapiro(data) if len(data) <= 5000 else (None, 0.001)
        return p_value > alpha if p_value is not None else False
    except:
        return False

# 3. BAND ÉS SCALE-SPECIFIKUS ELEMZÉS - NESTED FOR CIKLUS
print(f"\n2. Band és Scale-specifikus elemzés...")

# MOCA küszöb: -2 (mint az eredeti kódban)
moca_threshold = -2

# Minden eredmény tárolása
all_band_scale_results = []

# NESTED FOR CIKLUS - MINDEN SCALE ÉS BAND KOMBINÁCIÓRA
combination_counter = 0

for scale_idx, current_scale in enumerate(unique_scales):
    for band_idx, current_band in enumerate(unique_bands):
        combination_counter += 1
        
        print(f"\n{'='*70}")
        print(f"=== KOMBINÁCIÓ {combination_counter}/{total_combinations}: {current_scale} + {current_band} ===")
        print(f"{'='*70}")
        
        try:
            # Csak az aktuális scale ÉS band adatai
            combo_biomarker_df = biomarker_df[
                (biomarker_df['Scale'] == current_scale) & 
                (biomarker_df['Band'] == current_band)
            ].copy()
            
            print(f"Scale {current_scale} + Band {current_band} sorok: {len(combo_biomarker_df)}")
            
            # Merge MOCA adatokkal
            combo_merged_df = combo_biomarker_df.merge(moca_selected, left_on='PatientID', right_on='Identifier', how='inner')
            print(f"Merge után: {len(combo_merged_df)} sor")
            
            # Valid MOCA adatok
            valid_moca_mask = combo_merged_df['DeltaMOCA1'].notna()
            combo_final_df = combo_merged_df[valid_moca_mask].copy()
            print(f"Valid MOCA adatok: {len(combo_final_df)} sor (ez most tényleg unique betegek!)")
            
            if len(combo_final_df) < 10:
                print(f"SKIP {current_scale}+{current_band}: Túl kevés valid adat!")
                continue
            
            # MOCA csoportosítás erre a scale+band kombinációra
            good_moca_mask = combo_final_df['DeltaMOCA1'] >= moca_threshold
            bad_moca_mask = combo_final_df['DeltaMOCA1'] < moca_threshold
            
            good_moca_group = combo_final_df[good_moca_mask]
            bad_moca_group = combo_final_df[bad_moca_mask]
            
            print(f"Jó MOCA csoport (≥{moca_threshold}): {len(good_moca_group)} páciens")
            print(f"Rossz MOCA csoport (<{moca_threshold}): {len(bad_moca_group)} páciens")
            
            if len(good_moca_group) < 3 or len(bad_moca_group) < 3:
                print(f"SKIP {current_scale}+{current_band}: Túl kevés adat a csoportokban!")
                continue
            
            # Metrikák azonosítása
            exclude_columns = ['PatientID', 'Band', 'Scale', 'Identifier', 'DeltaMOCA1']
            metric_columns = [col for col in combo_final_df.columns if col not in exclude_columns]
            print(f"Elemzendő metrikák: {len(metric_columns)}")
            
            # Metrika típusok eloszlása
            metric_types = {}
            for col in metric_columns:
                metric_types[col] = detect_metric_type(col)
            
            # Statisztikai elemzés erre a scale+band kombinációra
            print(f"\nStatisztikai elemzés {len(metric_columns)} metrikára...")
            
            combo_results = []
            
            for idx, metric_name in enumerate(metric_columns):
                if idx < 5 or idx % 20 == 0:  # Csak első 5-öt és minden 20.-at printeljük (túl sok output)
                    print(f"\n--- {idx+1}/{len(metric_columns)}: {metric_name} ---")
                
                # Metrika típusa
                metric_type = metric_types[metric_name]
                if idx < 5 or idx % 20 == 0:
                    print(f"  Típus: {metric_type}")
                
                # Adatok kinyerése
                good_metric_data = good_moca_group[metric_name]
                bad_metric_data = bad_moca_group[metric_name]
                
                # Valid adatok száma
                good_valid = good_metric_data.notna().sum()
                bad_valid = bad_metric_data.notna().sum()
                if idx < 5 or idx % 20 == 0:
                    print(f"  Valid adatok: jó MOCA n={good_valid}, rossz MOCA n={bad_valid}")
                
                if good_valid < 3 or bad_valid < 3:
                    if idx < 5 or idx % 20 == 0:
                        print(f"  SKIP: Túl kevés valid adat")
                    continue
                
                # Statisztikai teszt
                test_result = choose_statistical_test(metric_type, good_metric_data, bad_metric_data)
                
                if test_result[0] is None:
                    if idx < 5 or idx % 20 == 0:
                        print(f"  SKIP: Teszt sikertelen")
                    continue
                
                test_type, p_value, effect_size, normality_assumption, additional_info = test_result
                
                # Szignifikancia
                significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
                
                # Leíró statisztikák
                good_clean = good_metric_data.dropna()
                bad_clean = bad_metric_data.dropna()
                
                # Eredmény tárolása
                result = {
                    'Scale': current_scale,
                    'Band': current_band,
                    'Metric': metric_name,
                    'Metric_Type': metric_type,
                    'Test_Type': test_type,
                    'p_value': p_value,
                    'Effect_Size': effect_size,
                    'Significance': significance,
                    'Good_MOCA_N': len(good_clean),
                    'Good_MOCA_Mean': np.mean(good_clean),
                    'Good_MOCA_Std': np.std(good_clean, ddof=1),
                    'Good_MOCA_Median': np.median(good_clean),
                    'Bad_MOCA_N': len(bad_clean),
                    'Bad_MOCA_Mean': np.mean(bad_clean),
                    'Bad_MOCA_Std': np.std(bad_clean, ddof=1),
                    'Bad_MOCA_Median': np.median(bad_clean),
                    'Normality_Assumption': normality_assumption,
                    'Additional_Info': additional_info,
                    'Mean_Difference': np.mean(good_clean) - np.mean(bad_clean),
                    'Median_Difference': np.median(good_clean) - np.median(bad_clean)
                }
                
                combo_results.append(result)
                all_band_scale_results.append(result)
                
                if idx < 5 or idx % 20 == 0 or significance != "ns":
                    print(f"  {test_type}: p={p_value:.4f}{significance}, effect={effect_size:.3f}")
                    print(f"    Jó MOCA: μ={np.mean(good_clean):.4f}, σ={np.std(good_clean, ddof=1):.4f}")
                    print(f"    Rossz MOCA: μ={np.mean(bad_clean):.4f}, σ={np.std(bad_clean, ddof=1):.4f}")
            
            # Scale+Band specifikus eredmények mentése
            if len(combo_results) > 0:
                combo_results_df = pd.DataFrame(combo_results)
                
                # Összes eredmény erre a scale+band kombinációra
                combo_filename = f'biomarker_moca_analysis_{current_scale}_{current_band}.csv'
                combo_results_df.to_csv(combo_filename, index=False)
                print(f"\n{current_scale}+{current_band} eredmények: {combo_filename} ({len(combo_results_df)} metrika)")
                
                # Szignifikáns eredmények erre a scale+band kombinációra
                combo_significant_df = combo_results_df[combo_results_df['p_value'] < 0.05]
                if len(combo_significant_df) > 0:
                    combo_sig_filename = f'significant_biomarker_moca_{current_scale}_{current_band}.csv'
                    combo_significant_df.to_csv(combo_sig_filename, index=False)
                    print(f"{current_scale}+{current_band} szignifikáns: {combo_sig_filename} ({len(combo_significant_df)} metrika)")
                    
                    # Top 3 legszignifikánsabb
                    top_3 = combo_significant_df.nsmallest(3, 'p_value')
                    print(f"\nTop 3 legszignifikánsabb {current_scale}+{current_band}-on:")
                    for _, row in top_3.iterrows():
                        print(f"  {row['Metric']}: p={row['p_value']:.6f}{row['Significance']}, effect={row['Effect_Size']:.3f}")
                else:
                    print(f"{current_scale}+{current_band}: Nincs szignifikáns eredmény")
            
        except Exception as e:
            print(f"HIBA {current_scale}+{current_band} elemzésnél: {e}")
            continue

# 4. ÖSSZES EREDMÉNY ÖSSZEFOGLALÁSA
print(f"\n{'='*80}")
print(f"=== ÖSSZES BAND×SCALE KOMBINÁCIÓ EREDMÉNYEINEK ÖSSZEFOGLALÁSA ===")
print(f"{'='*80}")

if len(all_band_scale_results) > 0:
    # Összes eredmény DataFrame
    all_results_df = pd.DataFrame(all_band_scale_results)
    print(f"Összes elemzett metrika (minden band×scale): {len(all_results_df)}")
    
    # Band×Scale kombinációnkénti statisztikák
    print(f"\nBand×Scale kombinációnkénti összesítés:")
    for scale in unique_scales:
        for band in unique_bands:
            combo_results = all_results_df[
                (all_results_df['Scale'] == scale) & 
                (all_results_df['Band'] == band)
            ]
            if len(combo_results) > 0:
                combo_significant = combo_results[combo_results['p_value'] < 0.05]
                print(f"  {scale}+{band}: {len(combo_results)} metrika, {len(combo_significant)} szignifikáns ({len(combo_significant)/len(combo_results)*100:.1f}%)")
    
    # Scale-enkénti összesítés
    print(f"\nScale-enkénti összesítés (minden band összesen):")
    for scale in unique_scales:
        scale_results = all_results_df[all_results_df['Scale'] == scale]
        if len(scale_results) > 0:
            scale_significant = scale_results[scale_results['p_value'] < 0.05]
            print(f"  {scale}: {len(scale_results)} metrika, {len(scale_significant)} szignifikáns ({len(scale_significant)/len(scale_results)*100:.1f}%)")
    
    # Band-enkénti összesítés
    print(f"\nBand-enkénti összesítés (minden scale összesen):")
    for band in unique_bands:
        band_results = all_results_df[all_results_df['Band'] == band]
        if len(band_results) > 0:
            band_significant = band_results[band_results['p_value'] < 0.05]
            print(f"  {band}: {len(band_results)} metrika, {len(band_significant)} szignifikáns ({len(band_significant)/len(band_results)*100:.1f}%)")
    
    # Összes szignifikáns eredmény
    all_significant_df = all_results_df[all_results_df['p_value'] < 0.05]
    print(f"\nÖsszes szignifikáns eredmény: {len(all_significant_df)}/{len(all_results_df)} ({len(all_significant_df)/len(all_results_df)*100:.1f}%)")
    
    # Metrika típusonkénti bontás
    if len(all_significant_df) > 0:
        print(f"\nMetrika típusonkénti szignifikáns eredmények:")
        type_sig_counts = all_significant_df['Metric_Type'].value_counts()
        for mtype, count in type_sig_counts.items():
            total_type = all_results_df[all_results_df['Metric_Type'] == mtype].shape[0]
            print(f"  {mtype}: {count}/{total_type} ({count/total_type*100:.1f}%)")
    
    # Top 10 legjobb eredmény összesen
    if len(all_significant_df) > 0:
        print(f"\nTop 10 legjobb eredmény (minden kombináció):")
        top_10_overall = all_significant_df.nsmallest(10, 'p_value')
        for _, row in top_10_overall.iterrows():
            print(f"  {row['Scale']}+{row['Band']}: {row['Metric']} (p={row['p_value']:.6f}, effect={row['Effect_Size']:.3f})")
    
    # Legjobb eredmények kombinációnként
    if len(all_significant_df) > 0:
        print(f"\nLegjobb eredmény minden kombinációban:")
        for scale in unique_scales:
            for band in unique_bands:
                combo_sig = all_significant_df[
                    (all_significant_df['Scale'] == scale) & 
                    (all_significant_df['Band'] == band)
                ]
                if len(combo_sig) > 0:
                    best = combo_sig.loc[combo_sig['p_value'].idxmin()]
                    print(f"  {scale}+{band}: {best['Metric']} (p={best['p_value']:.6f}, effect={best['Effect_Size']:.3f})")
    
    # Fájlba mentés - ÖSSZESÍTETT
    print(f"\n=== FÁJLBA MENTÉS ===")
    
    # Összes eredmény (minden band×scale)
    all_results_df.to_csv('comprehensive_biomarker_moca_analysis_ALL_BAND_SCALE.csv', index=False)
    print(f"Összes eredmény: comprehensive_biomarker_moca_analysis_ALL_BAND_SCALE.csv ({len(all_results_df)} sor)")
    
    # Összes szignifikáns eredmény
    if len(all_significant_df) > 0:
        all_significant_df.to_csv('significant_biomarker_moca_results_ALL_BAND_SCALE.csv', index=False)
        print(f"Összes szignifikáns: significant_biomarker_moca_results_ALL_BAND_SCALE.csv ({len(all_significant_df)} sor)")
    
    # Scale-specifikus összesítések (minden band összevonva scale-enként)
    for scale in unique_scales:
        scale_all_bands = all_results_df[all_results_df['Scale'] == scale]
        if len(scale_all_bands) > 0:
            scale_filename = f'biomarker_moca_summary_{scale}_ALL_BANDS.csv'
            scale_all_bands.to_csv(scale_filename, index=False)
            print(f"Scale {scale} összesítés: {scale_filename} ({len(scale_all_bands)} sor)")
    
    # Band-specifikus összesítések (minden scale összevonva band-enként)
    for band in unique_bands:
        band_all_scales = all_results_df[all_results_df['Band'] == band]
        if len(band_all_scales) > 0:
            band_filename = f'biomarker_moca_summary_{band}_ALL_SCALES.csv'
            band_all_scales.to_csv(band_filename, index=False)
            print(f"Band {band} összesítés: {band_filename} ({len(band_all_scales)} sor)")
    
    # Típusonkénti bontás (minden band×scale-ről)
    for mtype in all_results_df['Metric_Type'].unique():
        type_df = all_results_df[all_results_df['Metric_Type'] == mtype]
        filename = f'biomarker_moca_{mtype}_results_ALL_BAND_SCALE.csv'
        type_df.to_csv(filename, index=False)
        print(f"{mtype} típus (összes band×scale): {filename} ({len(type_df)} sor)")
    
    # Erős hatású eredmények (effect size > 0.5)
    strong_effect_df = all_results_df[all_results_df['Effect_Size'] > 0.5]
    if len(strong_effect_df) > 0:
        strong_effect_df.to_csv('strong_effect_biomarker_moca_results_ALL_BAND_SCALE.csv', index=False)
        print(f"Erős hatás (összes band×scale): strong_effect_biomarker_moca_results_ALL_BAND_SCALE.csv ({len(strong_effect_df)} sor)")
    
    print(f"\n=== ELEMZÉS BEFEJEZVE ===")
    print(f"Elemzett kombinációk: {len(unique_scales)} scale × {len(unique_bands)} band = {total_combinations}")
    print(f"Összes elemzett metrika: {len(all_results_df)}")
    print(f"Összes szignifikáns: {len(all_significant_df)} ({len(all_significant_df)/len(all_results_df)*100:.1f}%)")
    
else:
    print("❌ Nincs eredmény - problémák az elemzéssel!")

print(f"\nKimeneti fájlok:")
print(f"Band×Scale specifikus fájlok:")
for scale in unique_scales:
    for band in unique_bands:
        print(f"  - biomarker_moca_analysis_{scale}_{band}.csv")
        print(f"  - significant_biomarker_moca_{scale}_{band}.csv (ha van szignifikáns)")

print(f"\nÖsszesített fájlok:")
print(f"  - comprehensive_biomarker_moca_analysis_ALL_BAND_SCALE.csv (összes)")
print(f"  - significant_biomarker_moca_results_ALL_BAND_SCALE.csv (szignifikáns)")

print(f"\nScale összesítések:")
for scale in unique_scales:
    print(f"  - biomarker_moca_summary_{scale}_ALL_BANDS.csv")

print(f"\nBand összesítések:")
for band in unique_bands:
    print(f"  - biomarker_moca_summary_{band}_ALL_SCALES.csv")

print(f"\nTípus szerinti:")
print(f"  - biomarker_moca_[type]_results_ALL_BAND_SCALE.csv")
print(f"  - strong_effect_biomarker_moca_results_ALL_BAND_SCALE.csv")
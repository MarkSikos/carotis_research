import pandas as pd
import numpy as np
import json

# CSV beolvasása
df = pd.read_csv('specific_indices_matrix.csv')

# PatientID oszlop eltávolítása
data = df.drop('PatientID', axis=1)

# Csak azokat a sorokat távolítjuk el, ahol MINDEN érték NaN
data_clean = data.dropna(how='all')

# Korrelációs mátrix számítása
correlation_matrix = data_clean.corr()

print("Korrelációs mátrix sikeresen kiszámítva!")
print(f"Mátrix mérete: {correlation_matrix.shape}")
print(f"Adatsorok száma: {len(data_clean)}")
print(f"Mutatók száma: {len(data_clean.columns)}")

# Alapstatisztikák
print(f"\nKorrelációs statisztikák:")
print(f"Átlag korreláció (diagonális nélkül): {correlation_matrix.values[np.triu_indices_from(correlation_matrix.values, k=1)].mean():.3f}")
print(f"Max korreláció (diagonális nélkül): {correlation_matrix.values[np.triu_indices_from(correlation_matrix.values, k=1)].max():.3f}")
print(f"Min korreláció (diagonális nélkül): {correlation_matrix.values[np.triu_indices_from(correlation_matrix.values, k=1)].min():.3f}")

# Korrelációs mátrix mentése CSV formátumban
correlation_matrix.to_csv('correlation_matrix.csv')
print("\n✅ Korrelációs mátrix elmentve: correlation_matrix.csv")

# Adatok előkészítése a vizualizációhoz
def prepare_for_visualization():
    """Készíti elő az adatokat a JavaScript vizualizációhoz"""
    
    # Korrelációs értékek listává alakítása
    n = len(correlation_matrix.columns)
    data_list = []
    x_labels = []
    y_labels = []
    
    for i, row_name in enumerate(correlation_matrix.index):
        for j, col_name in enumerate(correlation_matrix.columns):
            data_list.append(correlation_matrix.iloc[i, j])
            x_labels.append(col_name)
            y_labels.append(row_name)
    
    # JSON objektum létrehozása
    viz_data = {
        'data': data_list,
        'x': x_labels,
        'y': y_labels,
        'metrics': list(correlation_matrix.columns),
        'shape': [n, n]
    }
    
    return viz_data

# Vizualizációs adatok exportálása
viz_data = prepare_for_visualization()

# JSON fájlba mentés
with open('correlation_data.json', 'w', encoding='utf-8') as f:
    json.dump(viz_data, f, ensure_ascii=False, indent=2)

print("✅ Vizualizációs adatok elmentve: correlation_data.json")

# Top korrelációk keresése (diagonális elemek nélkül)
def find_top_correlations(correlation_matrix, n_top=10):
    """Megkeresi a legerősebb korrelációkat"""
    
    # Felső háromszög mátrix (diagonális nélkül)
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool), k=1)
    correlations = []
    
    for i in range(len(correlation_matrix.index)):
        for j in range(len(correlation_matrix.columns)):
            if mask[i, j] and not np.isnan(correlation_matrix.iloc[i, j]):
                correlations.append({
                    'metric1': correlation_matrix.index[i],
                    'metric2': correlation_matrix.columns[j],
                    'correlation': correlation_matrix.iloc[i, j]
                })
    
    # Rendezés abszolút érték szerint
    correlations.sort(key=lambda x: abs(x['correlation']), reverse=True)
    
    return correlations[:n_top]

# Top korrelációk megjelenítése
print(f"\n🔝 Top 10 legerősebb korreláció:")
top_correlations = find_top_correlations(correlation_matrix, 10)

for i, corr in enumerate(top_correlations, 1):
    print(f"{i:2d}. {corr['metric1'][:30]:<30} ↔ {corr['metric2'][:30]:<30} | r = {corr['correlation']:6.3f}")

# Kategóriák szerinti csoportosítás
def categorize_metrics(metrics):
    """Kategorizálja a mutatókat"""
    categories = {
        'Cardiac': [],
        'Endothelial': [],
        'Neurogenic': [],
        'Myogenic': [],
        'Respiratory': [],
        'VLF': [],
        'LF': [],
        'Global': []
    }
    
    for metric in metrics:
        if 'Cardiac' in metric:
            categories['Cardiac'].append(metric)
        elif 'Endothelial' in metric:
            categories['Endothelial'].append(metric)
        elif 'Neurogenic' in metric:
            categories['Neurogenic'].append(metric)
        elif 'Myogenic' in metric:
            categories['Myogenic'].append(metric)
        elif 'Respiratory' in metric:
            categories['Respiratory'].append(metric)
        elif 'VLF' in metric:
            categories['VLF'].append(metric)
        elif 'LF' in metric:
            categories['LF'].append(metric)
        else:
            categories['Global'].append(metric)
    
    return categories

# Kategória statisztikák
categories = categorize_metrics(correlation_matrix.columns)
print(f"\n📊 Mutatók kategóriák szerint:")
for cat, metrics in categories.items():
    if metrics:
        print(f"{cat:12s}: {len(metrics):3d} mutató")

# Hiányzó adatok elemzése
print(f"\n❌ Hiányzó adatok elemzése:")
missing_data = data_clean.isnull().sum()
missing_data = missing_data[missing_data > 0].sort_values(ascending=False)

if len(missing_data) > 0:
    print(f"Hiányzó adatokat tartalmazó mutatók száma: {len(missing_data)}")
    print("Top 5 legtöbb hiányzó adattal:")
    for metric, missing_count in missing_data.head().items():
        missing_pct = (missing_count / len(data_clean)) * 100
        print(f"  {metric[:50]:<50} | {missing_count:3d} hiányzó ({missing_pct:5.1f}%)")
else:
    print("Nincsenek hiányzó adatok a megtisztított adathalmazban.")

print(f"\n✨ Minden fájl sikeresen létrehozva!")
print(f"📁 Fájlok:")
print(f"   - correlation_matrix.csv      (teljes korrelációs mátrix)")
print(f"   - correlation_data.json       (vizualizációs adatok)")
print(f"\n🌐 A HTML vizualizáció megnyitásához töltsd be a correlation_data.json fájlt!")

# Opcionális: Plotly heatmap létrehozása
try:
    import plotly.graph_objects as go
    import plotly.express as px
    
    # Egyszerű heatmap
    fig = go.Figure(data=go.Heatmap(
        z=correlation_matrix.values,
        x=correlation_matrix.columns,
        y=correlation_matrix.index,
        colorscale='RdBu',
        zmid=0,
        zmin=-1,
        zmax=1,
        colorbar=dict(title="Korreláció")
    ))
    
    fig.update_layout(
        title="Cerebrális Autoregulációs Metrikák - Korrelációs Mátrix",
        xaxis_title="Mutatók",
        yaxis_title="Mutatók",
        width=1200,
        height=1200
    )
    
    # HTML fájlba mentés
    fig.write_html("correlation_heatmap_plotly.html")
    print(f"   - correlation_heatmap_plotly.html (Plotly heatmap)")
    
except ImportError:
    print("\n📝 Megjegyzés: Plotly nincs telepítve. Ha szeretnéd a Plotly verziót:")
    print("   pip install plotly")
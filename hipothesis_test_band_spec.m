%% SÁV-SPECIFIKUS METRIKÁK ÉS MOCA ANALÍZIS
% Minden sáv minden metrikájára Down vs Same összehasonlítás

clear; clc;

fprintf('=== SÁV-SPECIFIKUS MOCA ANALÍZIS ===\n');

%% 1. Eredeti adatok (MOCA-val) betöltése
full_data = readtable('df_unnormalized.csv');

% MOCA csoportok kialakítása
down_patients = full_data.Identifier(full_data.DeltaMOCA1 <= -2);  % Kognitív romlás
same_patients = full_data.Identifier(full_data.DeltaMOCA1 > -2);   % Stabil/javuló

down_patients = unique(down_patients(~isnan(down_patients)));
same_patients = unique(same_patients(~isnan(same_patients)));

fprintf('MOCA csoportok:\n');
fprintf('  DOWN (≤-3): %d páciens\n', length(down_patients));
fprintf('  SAME (>-3): %d páciens\n', length(same_patients));

%% 2. Sáv-specifikus eredmények betöltése
band_names = {'Endothelial', 'Neurogenic', 'Myogenic', 'Respiratory', 'Cardiac'};
band_data = struct();

for i = 1:length(band_names)
    band_name = band_names{i};
    filename = sprintf('band_%s_patient_metrics.csv', lower(band_name));
    
    if exist(filename, 'file')
        band_data.(band_name) = readtable(filename);
        fprintf('  %s: %d sor betöltve (%s)\n', band_name, height(band_data.(band_name)), filename);
    else
        fprintf('  %s: FÁJL NEM TALÁLHATÓ (%s)\n', band_name, filename);
        band_data.(band_name) = table(); % Üres táblázat
    end
end

%% 3. HIPOTÉZIS TESZTELÉS - Minden sáv minden metrikájára
fprintf('\n=== HIPOTÉZIS TESZTELÉS ===\n');

% Összesített eredmény táblázat
all_test_results = table();

for i = 1:length(band_names)
    band_name = band_names{i};
    current_data = band_data.(band_name);
    
    if isempty(current_data)
        fprintf('\n--- %s SÁV: NINCS ADAT ---\n', band_name);
        continue;
    end
    
    fprintf('\n--- %s SÁV HIPOTÉZIS TESZTEK ---\n', band_name);
    
    % Dinamikusan felismerni a metrikákat (PatientID és Band kihagyása)
    all_columns = current_data.Properties.VariableNames;
    metric_columns = all_columns(~ismember(all_columns, {'PatientID', 'Band'}));
    
    fprintf('  Metrikák: %s\n', strjoin(metric_columns, ', '));
    
    % Minden metrikára tesztelés
    for m = 1:length(metric_columns)
        metric_name = metric_columns{m};
        metric_data = current_data.(metric_name);
        
        % Csoportokra bontás
        patient_ids = current_data.PatientID;
        down_indices = ismember(patient_ids, down_patients);
        same_indices = ismember(patient_ids, same_patients);
        
        down_data = metric_data(down_indices);
        same_data = metric_data(same_indices);
        
        % Csak valid adatok
        down_valid = down_data(~isnan(down_data) & isfinite(down_data));
        same_valid = same_data(~isnan(same_data) & isfinite(same_data));
        
        % Elégséges adat ellenőrzés
        if length(down_valid) >= 3 && length(same_valid) >= 3
            
            % Normalitás ellenőrzés mindkét csoportban
            down_normal = false;
            same_normal = false;
            
            if length(down_valid) >= 3 && length(down_valid) <= 50
                try
                    [~, down_p_norm] = swtest(down_valid);
                    down_normal = (down_p_norm >= 0.05);
                catch
                    down_normal = false; % Ha swtest hibázik
                end
            end
            
            if length(same_valid) >= 3 && length(same_valid) <= 50
                try
                    [~, same_p_norm] = swtest(same_valid);
                    same_normal = (same_p_norm >= 0.05);
                catch
                    same_normal = false;
                end
            end
            
            % Teszt kiválasztása és végrehajtása
            try
                if down_normal && same_normal
                    % Mindkét csoport normális → t-test
                    [~, p_value, ~, stats] = ttest2(down_valid, same_valid);
                    test_type = 't-test';
                    
                    % Effect size (Cohen's d)
                    pooled_std = sqrt(((length(down_valid)-1)*var(down_valid) + (length(same_valid)-1)*var(same_valid)) / ...
                                    (length(down_valid) + length(same_valid) - 2));
                    if pooled_std > 0
                        cohens_d = (mean(down_valid) - mean(same_valid)) / pooled_std;
                        effect_size = abs(cohens_d);
                    else
                        effect_size = 0;
                    end
                    
                    test_stat = stats.tstat;
                    
                else
                    % Legalább egy nem normális → Wilcoxon ranksum
                    [p_value, ~, stats] = ranksum(down_valid, same_valid);
                    test_type = 'Wilcoxon';
                    
                    % Effect size (r = Z/sqrt(N))
                    z_score = stats.zval;
                    total_n = length(down_valid) + length(same_valid);
                    effect_size = abs(z_score) / sqrt(total_n);
                    test_stat = z_score;
                end
                
                % Szignifikancia szint
                significance = '';
                if p_value < 0.001
                    significance = '***';
                elseif p_value < 0.01
                    significance = '**';
                elseif p_value < 0.05
                    significance = '*';
                else
                    significance = 'ns';
                end
                
                % Hatás méret kategorizálás
                if strcmp(test_type, 't-test')
                    % Cohen's d kategóriák
                    if effect_size < 0.2
                        effect_cat = 'kis';
                    elseif effect_size < 0.5
                        effect_cat = 'közepes';
                    elseif effect_size < 0.8
                        effect_cat = 'nagy';
                    else
                        effect_cat = 'nagyon nagy';
                    end
                else
                    % r effect size kategóriák
                    if effect_size < 0.1
                        effect_cat = 'kis';
                    elseif effect_size < 0.3
                        effect_cat = 'közepes';
                    elseif effect_size < 0.5
                        effect_cat = 'nagy';
                    else
                        effect_cat = 'nagyon nagy';
                    end
                end
                
                % Kompakt kimenet
                fprintf('  %s: %s, p=%.4f%s, effect=%.3f(%s), DOWN_μ=%.3f, SAME_μ=%.3f\n', ...
                       metric_name, test_type, p_value, significance, effect_size, effect_cat, ...
                       mean(down_valid), mean(same_valid));
                
                % Eredmény táblázatba mentés
                new_test_row = table({band_name}, {metric_name}, {test_type}, ...
                                    p_value, test_stat, effect_size, {effect_cat}, {significance}, ...
                                    length(down_valid), mean(down_valid), std(down_valid), ...
                                    length(same_valid), mean(same_valid), std(same_valid), ...
                                    'VariableNames', {'Band', 'Metric', 'Test_Type', 'p_value', 'Test_Stat', ...
                                    'Effect_Size', 'Effect_Category', 'Significance', ...
                                    'DOWN_N', 'DOWN_Mean', 'DOWN_Std', 'SAME_N', 'SAME_Mean', 'SAME_Std'});
                all_test_results = [all_test_results; new_test_row];
                
            catch ME
                fprintf('  %s: TESZT HIBA - %s\n', metric_name, ME.message);
            end
            
        else
            fprintf('  %s: Elégtelenül adat (DOWN n=%d, SAME n=%d)\n', metric_name, length(down_valid), length(same_valid));
        end
    end
end

%% 4. SZIGNIFIKÁNS EREDMÉNYEK ÖSSZEGZÉSE
fprintf('\n=== SZIGNIFIKÁNS EREDMÉNYEK ÖSSZEGZÉSE ===\n');

% Szignifikáns eredmények (p < 0.05)
significant_results = all_test_results(all_test_results.p_value < 0.05, :);

if height(significant_results) > 0
    fprintf('Szignifikáns különbségek: %d/%d teszt\n', height(significant_results), height(all_test_results));
    
    % Sávonkénti bontás
    for i = 1:length(band_names)
        band_name = band_names{i};
        band_significant = significant_results(strcmp(significant_results.Band, band_name), :);
        
        if height(band_significant) > 0
            fprintf('\n--- %s SÁV SZIGNIFIKÁNS EREDMÉNYEK ---\n', band_name);
            
            for j = 1:height(band_significant)
                row = band_significant(j, :);
                direction = '';
                if row.DOWN_Mean > row.SAME_Mean
                    direction = 'DOWN > SAME';
                else
                    direction = 'SAME > DOWN';
                end
                
                fprintf('  %s: %s, p=%.4f%s, %s effect, %s\n', ...
                       row.Metric{1}, row.Test_Type{1}, ...
                       row.p_value, row.Significance{1}, row.Effect_Category{1}, direction);
            end
        end
    end
    
    % Nagy hatás méret eredmények
    large_effect = significant_results(strcmp(significant_results.Effect_Category, 'nagy') | ...
                                     strcmp(significant_results.Effect_Category, 'nagyon nagy'), :);
    if height(large_effect) > 0
        fprintf('\n=== NAGY HATÁS MÉRETŰ KÜLÖNBSÉGEK ===\n');
        for i = 1:height(large_effect)
            row = large_effect(i, :);
            fprintf('  %s-%s: effect=%.3f (%s), p=%.4f%s\n', ...
                   row.Band{1}, row.Metric{1}, row.Effect_Size, row.Effect_Category{1}, ...
                   row.p_value, row.Significance{1});
        end
    end
    
else
    fprintf('Nincs szignifikáns különbség egyetlen metrikában sem (p ≥ 0.05)\n');
end

%% 5. STATISZTIKAI ÖSSZEFOGLALÓK
fprintf('\n=== STATISZTIKAI ÖSSZEFOGLALÓK ===\n');

% Teszt típusok összegzése
test_type_summary = groupcounts(all_test_results, 'Test_Type');
fprintf('Teszt típusok:\n');
disp(test_type_summary);

% Sávonkénti összefoglaló
band_summary = groupcounts(all_test_results, 'Band');
fprintf('\nSávonkénti tesztek száma:\n');
disp(band_summary);

% Szignifikancia szintek
sig_summary = groupcounts(all_test_results, 'Significance');
fprintf('\nSzignifikancia szintek:\n');
disp(sig_summary);

%% 6. EREDMÉNYEK EXPORTÁLÁSA
fprintf('\n=== EREDMÉNYEK EXPORTÁLÁSA ===\n');

% Összes teszt eredmény
writetable(all_test_results, 'band_specific_moca_tests.csv');
fprintf('Összes teszt eredmény: band_specific_moca_tests.csv (%d teszt)\n', height(all_test_results));

% Csak szignifikáns eredmények
if height(significant_results) > 0
    writetable(significant_results, 'band_specific_moca_significant.csv');
    fprintf('Szignifikáns eredmények: band_specific_moca_significant.csv (%d eredmény)\n', height(significant_results));
end

% Nagy hatás méretű eredmények
if height(large_effect) > 0
    writetable(large_effect, 'band_specific_moca_large_effects.csv');
    fprintf('Nagy hatás méretű eredmények: band_specific_moca_large_effects.csv (%d eredmény)\n', height(large_effect));
end

%% 7. KLINIKAI ÉRTELMEZÉS
fprintf('\n=== KLINIKAI ÉRTELMEZÉS ===\n');

if height(significant_results) > 0
    % Érintett sávok
    affected_bands = unique(significant_results.Band);
    fprintf('Érintett frekvencia sávok: %s\n', strjoin(affected_bands, ', '));
    
    % Érintett metrikák
    affected_metrics = unique(significant_results.Metric);
    fprintf('Érintett metrikák: %s\n', strjoin(affected_metrics, ', '));
    
    % Irányultság elemzés
    down_higher = sum(significant_results.DOWN_Mean > significant_results.SAME_Mean);
    same_higher = sum(significant_results.SAME_Mean > significant_results.DOWN_Mean);
    
    fprintf('\nIrányultság:\n');
    fprintf('  Kognitív romlás csoport magasabb: %d metrika\n', down_higher);
    fprintf('  Stabil csoport magasabb: %d metrika\n', same_higher);
    
    if down_higher > same_higher
        fprintf('  → Túlnyomórészt: Kognitív romlás magasabb értékekkel jár együtt\n');
    elseif same_higher > down_higher
        fprintf('  → Túlnyomórészt: Stabil kognitív állapot magasabb értékekkel jár együtt\n');
    else
        fprintf('  → Vegyes irányultság - részletes sávonkénti elemzés szükséges\n');
    end
    
    % Sávonkénti irányultság
    fprintf('\nSávonkénti irányultság:\n');
    for i = 1:length(affected_bands)
        band_name = affected_bands{i};
        band_results = significant_results(strcmp(significant_results.Band, band_name), :);
        
        band_down_higher = sum(band_results.DOWN_Mean > band_results.SAME_Mean);
        band_same_higher = sum(band_results.SAME_Mean > band_results.DOWN_Mean);
        
        fprintf('  %s: DOWN>SAME: %d, SAME>DOWN: %d\n', band_name, band_down_higher, band_same_higher);
    end
    
else
    fprintf('Nincs szignifikáns különbség → További elemzés vagy nagyobb minta szükséges\n');
end

fprintf('\n=== ANALÍZIS BEFEJEZVE ===\n');
fprintf('Kimeneti fájlok:\n');
fprintf('  - band_specific_moca_tests.csv (összes teszt)\n');
fprintf('  - band_specific_moca_significant.csv (szignifikáns)\n');
fprintf('  - band_specific_moca_large_effects.csv (nagy hatás)\n');
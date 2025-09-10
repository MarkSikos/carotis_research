%% WAVELET EREDMÉNYEK ÉS MOCA ANALÍZIS

% 1. Wavelet eredmények betöltése
wavelet_results = readtable('patient_wavelet_results.csv');
fprintf('Wavelet eredmények betöltve: %d sor\n', height(wavelet_results));

% 2. Eredeti adatok (MOCA-val) betöltése
full_data = readtable('df_unnormalized.csv');

% 3. MOCA csoportok kialakítása
down_patients = full_data.Identifier(full_data.DeltaMOCA1 <= -3);  % Kognitív romlás
same_patients = full_data.Identifier(full_data.DeltaMOCA1 > -3);   % Stabil/javuló

down_patients = unique(down_patients(~isnan(down_patients)));
same_patients = unique(same_patients(~isnan(same_patients)));

% 4. Wavelet eredmények szűrése - csak sikeres elemzések
wavelet_results = wavelet_results(~isnan(wavelet_results.Coherence), :);
% 5. Közös betegek azonosítása
wavelet_patients = unique(wavelet_results.PatientID);


%% HIPOTÉZIS TESZTELÉS - DOWN vs SAME CSOPORTOK ÖSSZEHASONLÍTÁSA

fprintf('\n=== HIPOTÉZIS TESZTELÉS ===\n');

% Frekvencia sávok és paraméterek
bands = unique(wavelet_results.Band);
params = {'Coherence', 'PLV', 'CMP_deg', 'COx'};

% Eredmény táblázat a tesztekhez
test_results = table();

for b = 1:length(bands)
    band_name = bands{b};
    band_data = wavelet_results(strcmp(wavelet_results.Band, band_name), :);
    
    fprintf('\n--- %s SÁV HIPOTÉZIS TESZTEK ---\n', band_name);
    
    for p = 1:length(params)
        param_name = params{p};
        param_data = band_data.(param_name);
        
        % Csoportokra bontás
        down_data = param_data(ismember(band_data.PatientID, down_patients));
        same_data = param_data(ismember(band_data.PatientID, same_patients));
        
        % Csak valid adatok
        down_valid = down_data(~isnan(down_data));
        same_valid = same_data(~isnan(same_data));
        
        if length(down_valid) >= 3 && length(same_valid) >= 3
            
            % Normalitás ellenőrzés mindkét csoportban
            down_normal = false;
            same_normal = false;
            
            if length(down_valid) >= 3 && length(down_valid) <= 50
                [~, down_p_norm] = swtest(down_valid);
                down_normal = (down_p_norm >= 0.05);
            end
            
            if length(same_valid) >= 3 && length(same_valid) <= 50
                [~, same_p_norm] = swtest(same_valid);
                same_normal = (same_p_norm >= 0.05);
            end
            
            % Teszt kiválasztása
            if down_normal && same_normal
                % Mindkét csoport normális → t-test
                [~, p_value, ~, stats] = ttest2(down_valid, same_valid);
                test_type = 't-test';
                
                % Effect size (Cohen's d)
                pooled_std = sqrt(((length(down_valid)-1)*var(down_valid) + (length(same_valid)-1)*var(same_valid)) / ...
                                (length(down_valid) + length(same_valid) - 2));
                cohens_d = (mean(down_valid) - mean(same_valid)) / pooled_std;
                effect_size = abs(cohens_d);
                
                % t-statisztika
                test_stat = stats.tstat;
                
            else
                % Legalább egy nem normális → Wilcoxon ranksum (Mann-Whitney U)
                [p_value, ~, stats] = ranksum(down_valid, same_valid);
                test_type = 'Wilcoxon';
                
                % Effect size (r = Z/sqrt(N))
                z_score = stats.zval;
                total_n = length(down_valid) + length(same_valid);
                effect_size = abs(z_score) / sqrt(total_n);
                
                % Z-statisztika
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
                % Cohen's d: 0.2=small, 0.5=medium, 0.8=large
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
                % r effect size: 0.1=small, 0.3=medium, 0.5=large
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
            fprintf('%s: %s, p=%.4f%s, effect=%.3f(%s), DOWN_μ=%.3f, SAME_μ=%.3f\n', ...
                   param_name, test_type, p_value, significance, effect_size, effect_cat, ...
                   mean(down_valid), mean(same_valid));
            
            % Eredmény táblázatba mentés
            new_test_row = table({band_name}, {param_name}, {test_type}, ...
                                p_value, test_stat, effect_size, {effect_cat}, {significance}, ...
                                length(down_valid), mean(down_valid), std(down_valid), ...
                                length(same_valid), mean(same_valid), std(same_valid), ...
                                'VariableNames', {'Band', 'Parameter', 'Test_Type', 'p_value', 'Test_Stat', ...
                                'Effect_Size', 'Effect_Category', 'Significance', ...
                                'DOWN_N', 'DOWN_Mean', 'DOWN_Std', 'SAME_N', 'SAME_Mean', 'SAME_Std'});
            test_results = [test_results; new_test_row];
            
        else
            fprintf('%s: Elégtelenül adat (DOWN n=%d, SAME n=%d)\n', param_name, length(down_valid), length(same_valid));
        end
    end
end

%% SZIGNIFIKÁNS EREDMÉNYEK ÖSSZEGZÉSE
fprintf('\n=== SZIGNIFIKÁNS EREDMÉNYEK ===\n');

% Szignifikáns eredmények (p < 0.05)
significant_results = test_results(test_results.p_value < 0.05, :);

if height(significant_results) > 0
    fprintf('Szignifikáns különbségek: %d/%d teszt\n', height(significant_results), height(test_results));
    
    for i = 1:height(significant_results)
        row = significant_results(i, :);
        direction = '';
        if row.DOWN_Mean > row.SAME_Mean
            direction = 'DOWN > SAME';
        else
            direction = 'SAME > DOWN';
        end
        
        fprintf('  %s-%s: %s, p=%.4f%s, %s effect, %s\n', ...
               row.Band{1}, row.Parameter{1}, row.Test_Type{1}, ...
               row.p_value, row.Significance{1}, row.Effect_Category{1}, direction);
    end
    
    % Nagy hatás méret eredmények
    large_effect = significant_results(strcmp(significant_results.Effect_Category, 'nagy') | ...
                                     strcmp(significant_results.Effect_Category, 'nagyon nagy'), :);
    if height(large_effect) > 0
        fprintf('\nNagy hatás méretű különbségek: %d\n', height(large_effect));
        for i = 1:height(large_effect)
            row = large_effect(i, :);
            fprintf('  %s-%s: effect=%.3f (%s)\n', row.Band{1}, row.Parameter{1}, ...
                   row.Effect_Size, row.Effect_Category{1});
        end
    end
    
else
    fprintf('Nincs szignifikáns különbség egyetlen paraméterben sem (p ≥ 0.05)\n');
end

%% TESZT TÍPUSOK ÖSSZEGZÉSE
fprintf('\n=== TESZT TÍPUSOK ===\n');
test_type_summary = groupcounts(test_results, 'Test_Type');
disp(test_type_summary);

%% EXPORT EREDMÉNYEK
writetable(test_results, 'moca_hypothesis_tests.csv');
fprintf('\nHipotézis teszt eredmények mentve: moca_hypothesis_tests.csv\n');

% Csak szignifikáns eredmények külön fájlba
if height(significant_results) > 0
    writetable(significant_results, 'moca_significant_results.csv');
    fprintf('Szignifikáns eredmények mentve: moca_significant_results.csv\n');
end

%% KLINIKAI ÉRTELMEZÉS
fprintf('\n=== KLINIKAI ÉRTELMEZÉS ===\n');

if height(significant_results) > 0
    % Mely sávokban vannak különbségek
    affected_bands = unique(significant_results.Band);
    fprintf('Érintett frekvencia sávok: %s\n', strjoin(affected_bands, ', '));
    
    % Mely paraméterekben vannak különbségek  
    affected_params = unique(significant_results.Parameter);
    fprintf('Érintett paraméterek: %s\n', strjoin(affected_params, ', '));
    
    % Irányultság
    down_higher = sum(significant_results.DOWN_Mean > significant_results.SAME_Mean);
    same_higher = sum(significant_results.SAME_Mean > significant_results.DOWN_Mean);
    
    fprintf('Kognitív romlás csoport magasabb: %d esetben\n', down_higher);
    fprintf('Stabil csoport magasabb: %d esetben\n', same_higher);
    
    if down_higher > same_higher
        fprintf('→ Kognitív romlás magasabb neurovascular coupling értékekkel jár együtt\n');
    elseif same_higher > down_higher
        fprintf('→ Stabil kognitív állapot magasabb neurovascular coupling értékekkel jár együtt\n');
    else
        fprintf('→ Vegyes irányultság - részletes elemzés szükséges\n');
    end
    
else
    fprintf('Nincs szignifikáns különbség → További elemzés vagy nagyobb minta szükséges\n');
end

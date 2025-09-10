clear; clc;

%% SÁV-SPECIFIKUS FEJLETT SZÖGELOSZLÁSOS ELEMZÉS
% A működő kód sávokra bontva + javított datetime kezelés

fprintf('=== SÁV-SPECIFIKUS FEJLETT SZÖGELOSZLÁSOS ELEMZÉS ===\n');

%% 1. PARAMÉTEREK ÉS FREKVENCIASÁVOK (RELAXÁLT)
data_file = 'data/df_unnormalized.csv';
target_oppart = 4;  % Clamp fázis
min_valid_ratio = 0.50;  % 🔧 RELAXÁLT: 70% -> 50%

% Deriváltfüggvény paraméterek (ugyanazok mint a működő verzióban)
derivative_methods = {'gradient', 'local_window', 'savitzky_golay', 'wavelet_based'};
window_sizes = [3, 5, 10, 20];  % Multi-scale elemzéshez
kde_bandwidth = 5;              % Kernel density estimation sávszélesség

% 🔧 RELAXÁLT FIZIOLÓGIAI HATÁROK
physio_bounds = [-20, 300; -5, 110; 0, 250]; % Bővített határok

%% 🎯 FREKVENCIASÁVOK DEFINÍCIÓJA (Csökkentett követelmények)
freq_bands = struct();

% Endothelial sáv - NO-mediált, metabolikus
freq_bands.Endothelial = struct( ...
    'range', [0.003, 0.02], ...
    'center_freq', 0.01, ...
    'target_fs', 0.1, ...
    'n_points', 200, ...           % 🔧 Csökkentett: 300->200
    'target_duration', 400, ...    % 🔧 Csökkentett: 800->400
    'mechanism', 'NO-mediated vasodilation, metabolic regulation', ...
    'clinical_relevance', 'Diabetes, endothelial dysfunction');

% Neurogenic sáv - szimpatikus, baroreflex
freq_bands.Neurogenic = struct( ...
    'range', [0.02, 0.06], ...
    'center_freq', 0.04, ...
    'target_fs', 0.2, ...
    'n_points', 250, ...           % 🔧 Csökkentett: 350->250
    'target_duration', 300, ...    % 🔧 Csökkentett: 600->300
    'mechanism', 'Sympathetic nervous system, baroreflex modulation', ...
    'clinical_relevance', 'Hypertension, autonomic dysfunction');

% Myogenic sáv - simaizom autoregulációs
freq_bands.Myogenic = struct( ...
    'range', [0.06, 0.15], ...
    'center_freq', 0.1, ...
    'target_fs', 0.5, ...
    'n_points', 300, ...           % 🔧 Csökkentett: 400->300
    'target_duration', 200, ...    % 🔧 Csökkentett: 400->200
    'mechanism', 'Smooth muscle autoregulation, Bayliss effect', ...
    'clinical_relevance', 'Arterial stiffness, aging, autoregulation');

% Aktív sávok
active_bands = {'Endothelial', 'Neurogenic', 'Myogenic'};
band_names = active_bands;

% Path hozzáadás
addpath(genpath('functions'));

%% 2. ADATOK BETÖLTÉSE ÉS SZŰRÉSE
fprintf('\n=== ADATOK BETÖLTÉSE ===\n');

data = readtable(data_file);
if ismember('oppart', data.Properties.VariableNames)
    phase_data = data(data.oppart == target_oppart, :);
elseif ismember('oper_phase', data.Properties.VariableNames)
    phase_data = data(data.oper_phase == target_oppart, :);
else
    fprintf('Elérhető oszlopok:\n');
    disp(data.Properties.VariableNames);
    error('Nem található oppart vagy oper_phase oszlop!');
end

patients = unique(phase_data.Identifier);
n_patients = length(patients);
fprintf('Eredeti páciensek száma: %d\n', n_patients);

%% 3. SÁV-SPECIFIKUS FELDOLGOZÁS MINDEN SÁVRA
all_band_results = struct();

for band_idx = 1:length(band_names)
    band_name = band_names{band_idx};
    band_params = freq_bands.(band_name);
    
    fprintf('\n==========================================\n');
    fprintf('=== %s SÁV ELEMZÉS ===\n', upper(band_name));
    fprintf('==========================================\n');
    fprintf('Frekvencia: %.3f-%.3f Hz (közép: %.3f Hz)\n', ...
            band_params.range(1), band_params.range(2), band_params.center_freq);
    fprintf('Mechanizmus: %s\n', band_params.mechanism);
    fprintf('Klinikai relevancia: %s\n', band_params.clinical_relevance);
    
    % Sáv-specifikus eredmények tárolása
    band_results = struct();
    band_results.band_name = band_name;
    band_results.band_params = band_params;
    band_results.patient_data = [];
    
    % 🔧 Sáv-specifikus gap tolerancia
    max_gap_points = round(60 * band_params.target_fs);  % 60 másodperc
    
    fprintf('\n--- %s sáv páciensenkénti feldolgozás ---\n', band_name);
    
    %% 3.1 PÁCIENSENKÉNTI FELDOLGOZÁS EZZEL A SÁVVAL
    for p = 1:n_patients
        patient_id = patients(p);
        patient_data = phase_data(phase_data.Identifier == patient_id, :);
        
        fprintf('Páciens %d feldolgozása (%s sáv)...\n', patient_id, band_name);
        
        %% RELAXÁLT ADATMINŐSÉG ELLENŐRZÉS
        SE_signal = patient_data.SE;
        rSO2_signal = patient_data.oper_side_oxig;
        MAP_signal = patient_data.MAP;
        rSO2_contra_signal = patient_data.other_side_oxig;
        
        se_valid_ratio = sum(~isnan(SE_signal)) / length(SE_signal);
        rso2_valid_ratio = sum(~isnan(rSO2_signal)) / length(rSO2_signal);
        map_valid_ratio = sum(~isnan(MAP_signal)) / length(MAP_signal);
        rso2_contra_valid_ratio = sum(~isnan(rSO2_contra_signal)) / length(rSO2_contra_signal);
        
        % 🔧 RELAXÁLT: 50% küszöb 70% helyett
        if se_valid_ratio < min_valid_ratio || rso2_valid_ratio < min_valid_ratio || ...
           map_valid_ratio < min_valid_ratio || rso2_contra_valid_ratio < min_valid_ratio
            fprintf('  ❌ Elégtelen adatminőség (SE:%.1f%%, rSO2:%.1f%%, MAP:%.1f%%)\n', ...
                    se_valid_ratio*100, rso2_valid_ratio*100, map_valid_ratio*100);
            continue;
        end
        
        try
            %% Gap filling (mint a működő verzióban)
            SE_cleaned = smart_gap_filling(SE_signal, max_gap_points);
            rSO2_cleaned = smart_gap_filling(rSO2_signal, max_gap_points);
            MAP_cleaned = smart_gap_filling(MAP_signal, max_gap_points);
            rSO2_contra_cleaned = smart_gap_filling(rSO2_contra_signal, max_gap_points);
            
            %% 🔥 JAVÍTOTT SÁVÁSPECIFIKUS INTERPOLÁCIÓ (Datetime-aware)
            patient_data_for_interp = patient_data;
            patient_data_for_interp.SE = SE_cleaned;
            patient_data_for_interp.oper_side_oxig = rSO2_cleaned;
            patient_data_for_interp.MAP = MAP_cleaned;
            patient_data_for_interp.other_side_oxig = rSO2_contra_cleaned;
            
            interpolated_data = interpolate_data_band_specific_fixed(patient_data_for_interp, band_params);
            
            SE_clean = fillmissing(interpolated_data.SE, 'nearest');
            rSO2_clean = fillmissing(interpolated_data.rSO2, 'nearest');
            MAP_clean = fillmissing(interpolated_data.MAP, 'nearest');
            rSO2_contra_clean = fillmissing(interpolated_data.rSO2_contra, 'nearest');
            
            %% RELAXÁLT FIZIOLÓGIAI ELLENŐRZÉS
            SE_physio_valid = SE_clean >= physio_bounds(1,1) & SE_clean <= physio_bounds(1,2);
            rSO2_physio_valid = rSO2_clean >= physio_bounds(2,1) & rSO2_clean <= physio_bounds(2,2);
            MAP_physio_valid = MAP_clean >= physio_bounds(3,1) & MAP_clean <= physio_bounds(3,2);
            rSO2_contra_physio_valid = rSO2_contra_clean >= physio_bounds(2,1) & rSO2_contra_clean <= physio_bounds(2,2);
            
            physio_valid_ratio = sum(SE_physio_valid & rSO2_physio_valid & MAP_physio_valid & rSO2_contra_physio_valid) / length(SE_clean);
            
            % 🔧 RELAXÁLT: 50% küszöb 70% helyett
            if physio_valid_ratio < min_valid_ratio
                fprintf('  ❌ Fiziológiai határok sértése (%.1f%% < %.0f%%)\n', physio_valid_ratio*100, min_valid_ratio*100);
                continue;
            end
            
            %% 🎯 SÁV-SPECIFIKUS SZŰRÉS ÉS FELDOLGOZÁS
            
            % Detrending (mint a működő verzióban)
            SE_detrended = detrend(SE_clean);
            rSO2_detrended = detrend(rSO2_clean);
            MAP_detrended = detrend(MAP_clean);                    % ← ÚJ
            rSO2_contra_detrended = detrend(rSO2_contra_clean); 
            
            % Sáv-specifikus bandpass szűrés HELYES actual_fs-sel
            [SE_filtered, rSO2_filtered] = apply_band_specific_filter_fixed(SE_detrended, rSO2_detrended, band_params, interpolated_data.actual_fs);
            [MAP_filtered, rSO2_contra_filtered] = apply_band_specific_filter_fixed(MAP_detrended, rSO2_contra_detrended, band_params, interpolated_data.actual_fs);
            %% 🔥 FEJLETT IRÁNYMENTI DERIVÁLTFÜGGVÉNY ELEMZÉS (SÁV-SPECIFIKUS)
            
            dt = interpolated_data.dt;  % 🔧 JAVÍTOTT: helyes dt használata
            time_vector = interpolated_data.time;  % 🔧 JAVÍTOTT: interpolált idővektor
            
            % Multi-scale derivált számítás (mint a működő verzióban)
            SE_derivatives = calculate_multiscale_derivatives(SE_filtered, dt, window_sizes);
            rSO2_derivatives = calculate_multiscale_derivatives(rSO2_filtered, dt, window_sizes);
            MAP_derivatives = calculate_multiscale_derivatives(MAP_filtered, dt, window_sizes);                      % ← ÚJ
            rSO2_contra_derivatives = calculate_multiscale_derivatives(rSO2_contra_filtered, dt, window_sizes); 

            % Kontinuus szögeloszlás minden skálán (mint a működő verzióban)
            SE_angle_analysis = struct();
            rSO2_angle_analysis = struct();
            MAP_angle_analysis = struct();                   % ← ÚJ
            rSO2_contra_angle_analysis = struct(); 
            
            for scale_idx = 1:length(window_sizes)
                scale_name = sprintf('scale_%d', window_sizes(scale_idx));
                
                % Szögek számítása ezen a skálán
                if isfield(SE_derivatives, ['scale_' num2str(window_sizes(scale_idx))])
                    SE_angles = atan2d(SE_derivatives.(['scale_' num2str(window_sizes(scale_idx))]), dt);
                    rSO2_angles = atan2d(rSO2_derivatives.(['scale_' num2str(window_sizes(scale_idx))]), dt);
                    MAP_angles = atan2d(MAP_derivatives.(['scale_' num2str(window_sizes(scale_idx))]), dt);                      % ← ÚJ
                    rSO2_contra_angles = atan2d(rSO2_contra_derivatives.(['scale_' num2str(window_sizes(scale_idx))]), dt);      % ← ÚJ
        
                    % KDE folytonos eloszlás (mint a működő verzióban)
                    [SE_kde, SE_kde_x] = calculate_kde(SE_angles, kde_bandwidth);
                    [rSO2_kde, rSO2_kde_x] = calculate_kde(rSO2_angles, kde_bandwidth);
                    [MAP_kde, MAP_kde_x] = calculate_kde(MAP_angles, kde_bandwidth);                      % ← ÚJ
                    [rSO2_contra_kde, rSO2_contra_kde_x] = calculate_kde(rSO2_contra_angles, kde_bandwidth);      % ← ÚJ    

                    % Fejlett statisztikai elemzés (mint a működő verzióban)
                    SE_analysis = analyze_continuous_angle_distribution(SE_angles, SE_kde, SE_kde_x);
                    rSO2_analysis = analyze_continuous_angle_distribution(rSO2_angles, rSO2_kde, rSO2_kde_x);
                    MAP_analysis = analyze_continuous_angle_distribution(MAP_angles, MAP_kde, MAP_kde_x);                      % ← ÚJ
                    rSO2_contra_analysis = analyze_continuous_angle_distribution(rSO2_contra_angles, rSO2_contra_kde, rSO2_contra_kde_x);      % ← ÚJ
                   
                    % Sáv-specifikus biomarkerek
                    SE_analysis.band_biomarkers = calculate_band_biomarkers(SE_analysis, band_name);
                    rSO2_analysis.band_biomarkers = calculate_band_biomarkers(rSO2_analysis, band_name);
                    MAP_analysis.band_biomarkers = calculate_band_biomarkers(MAP_analysis, band_name);                      % ← ÚJ
                    rSO2_contra_analysis.band_biomarkers = calculate_band_biomarkers(rSO2_contra_analysis, band_name);      % ← ÚJ
         
                    % Tárolás
                    SE_angle_analysis.(scale_name) = SE_analysis;
                    rSO2_angle_analysis.(scale_name) = rSO2_analysis;
                    MAP_angle_analysis.(scale_name) = MAP_analysis;                      % ← ÚJ
                    rSO2_contra_angle_analysis.(scale_name) = rSO2_contra_analysis;
                end
            end
            
            %% Wavelet-alapú frekvencia-függő szögeloszlás (mint a működő verzióban)
            if exist('cwt', 'file') == 2
                SE_wavelet_angles = calculate_wavelet_angle_distribution(SE_filtered, interpolated_data.actual_fs);
                rSO2_wavelet_angles = calculate_wavelet_angle_distribution(rSO2_filtered, interpolated_data.actual_fs);
                MAP_wavelet_angles = calculate_wavelet_angle_distribution(MAP_filtered, interpolated_data.actual_fs);                      % ← ÚJ
                rSO2_contra_wavelet_angles = calculate_wavelet_angle_distribution(rSO2_contra_filtered, interpolated_data.actual_fs);  
            else
                SE_wavelet_angles = [];
                rSO2_wavelet_angles = [];
                MAP_wavelet_angles = [];                      % ← ÚJ
                rSO2_contra_wavelet_angles = []; 
            end
            
            %% Keresztkorreláció szögeloszlások között (mint a működő verzióban)
            cross_angle_analysis = analyze_cross_angle_correlations(SE_angle_analysis, rSO2_angle_analysis);
            
            %% Sáv-specifikus spektrális teljesítmény
            [SE_band_power, rSO2_band_power] = calculate_band_power(SE_filtered, rSO2_filtered, band_params, interpolated_data.actual_fs);
            [MAP_band_power, rSO2_contra_band_power] = calculate_band_power(MAP_filtered, rSO2_contra_filtered, band_params, interpolated_data.actual_fs);  % ← ÚJ
            %% Eredmények tárolása (bővített a sáv infóval)
            %% Eredmények tárolása (bővített a sáv infóval) - BŐVÍTVE
            result = struct();
            result.PatientID = patient_id;
            result.band_name = band_name;
            result.SE_signal = SE_filtered;
            result.rSO2_signal = rSO2_filtered;
            result.MAP_signal = MAP_filtered;                          % ← ÚJ
            result.rSO2_contra_signal = rSO2_contra_filtered;          % ← ÚJ
            result.time_vector = time_vector;
            result.SE_derivatives = SE_derivatives;
            result.rSO2_derivatives = rSO2_derivatives;
            result.MAP_derivatives = MAP_derivatives;                  % ← ÚJ
            result.rSO2_contra_derivatives = rSO2_contra_derivatives;  % ← ÚJ
            result.SE_angle_analysis = SE_angle_analysis;
            result.rSO2_angle_analysis = rSO2_angle_analysis;
            result.MAP_angle_analysis = MAP_angle_analysis;                  % ← ÚJ
            result.rSO2_contra_angle_analysis = rSO2_contra_angle_analysis;  % ← ÚJ
            result.SE_wavelet_angles = SE_wavelet_angles;
            result.rSO2_wavelet_angles = rSO2_wavelet_angles;
            result.MAP_wavelet_angles = MAP_wavelet_angles;                  % ← ÚJ
            result.rSO2_contra_wavelet_angles = rSO2_contra_wavelet_angles;  % ← ÚJ
            result.cross_angle_analysis = cross_angle_analysis;
            result.SE_band_power = SE_band_power;
            result.rSO2_band_power = rSO2_band_power;
            result.MAP_band_power = MAP_band_power;                          % ← ÚJ
            result.rSO2_contra_band_power = rSO2_contra_band_power;          % ← ÚJ
            result.data_quality = physio_valid_ratio;
            result.signal_length = length(SE_filtered);
            result.actual_fs = interpolated_data.actual_fs;
            result.processing_timestamp = datetime('now');
            
            band_results.patient_data = [band_results.patient_data; result];
            
            fprintf('  ✅ Sikeres (%.1f%% valid, fs=%.3fHz, %d skála)\n', ...
                    physio_valid_ratio*100, interpolated_data.actual_fs, length(window_sizes));
            
        catch ME
            fprintf('  ❌ Feldolgozási hiba: %s\n', ME.message);
            continue;
        end
    end
    
    %% Sáv-specifikus összegzés
    n_successful = length(band_results.patient_data);
    fprintf('\n%s sáv összegzés:\n', band_name);
    fprintf('  Sikeres páciensek: %d/%d (%.1f%%)\n', n_successful, n_patients, n_successful/n_patients*100);
    
    if n_successful > 0
        % Sáv-specifikus populációs statisztikák
        band_results.population_stats = calculate_band_population_statistics(band_results.patient_data, band_name);
        
        % Sáv eredmények tárolása
        all_band_results.(band_name) = band_results;
        
        fprintf('  Populációs átlagok:\n');
        fprintf('    SE átlagszög: %.2f° ± %.2f°\n', ...
                band_results.population_stats.SE_mean_angle_mean, ...
                band_results.population_stats.SE_mean_angle_std);
        fprintf('    rSO2 átlagszög: %.2f° ± %.2f°\n', ...
                band_results.population_stats.rSO2_mean_angle_mean, ...
                band_results.population_stats.rSO2_mean_angle_std);
    else
        fprintf('  ⚠️ Nincs sikeres eredmény erre a sávra\n');
    end
end

%% 4. SÁVOK KÖZÖTTI ÖSSZEHASONLÍTÓ ELEMZÉS
successful_bands = fieldnames(all_band_results);
n_successful_bands = length(successful_bands);

fprintf('\n==========================================\n');
fprintf('=== SÁVOK KÖZÖTTI ÖSSZEHASONLÍTÁS ===\n');
fprintf('==========================================\n');
fprintf('Sikeresen elemzett sávok: %d/%d\n', n_successful_bands, length(band_names));

%% 5. KOMPREHENZÍV VIZUALIZÁCIÓ
if n_successful_bands > 0
    fprintf('\n=== SÁV-SPECIFIKUS VIZUALIZÁCIÓ ===\n');
    
    % MINDEN SIKERES SÁV VIZUALIZÁCIÓJA
    for i = 1:length(successful_bands)
        band_name = successful_bands{i};
        if ~isempty(all_band_results.(band_name).patient_data)
            example_patient = all_band_results.(band_name).patient_data(1);
            fprintf('Vizualizáció: %s sáv, Páciens %d\n', band_name, example_patient.PatientID);
            
            % Minden sáv külön figure-ben
            figure_title = sprintf('Sáv-%d-%s', i, band_name);
            set(0, 'CurrentFigure', figure('Name', figure_title, 'Position', [50+i*50, 50+i*50, 1400, 1000]));
            visualize_band_specific_analysis(example_patient, window_sizes, band_name);
            
            % Kis szünet a figure-k között
            pause(0.5);
        end
    end
    
    % ÖSSZEHASONLÍTÓ MULTI-SÁV VIZUALIZÁCIÓ
    if length(successful_bands) >= 2
        fprintf('\nÖsszehasonlító multi-sáv vizualizáció...\n');
        visualize_multi_band_comparison(all_band_results, successful_bands);
    end
end

fprintf('\n🎉 SÁV-SPECIFIKUS FEJLETT SZÖGELOSZLÁS ELEMZÉS BEFEJEZVE! 🎉\n');
fprintf('Sikeresen elemzett sávok: %s\n', strjoin(successful_bands, ', '));

%% ============================================================================
%% 🔥 JAVÍTOTT HELPER FÜGGVÉNYEK
%% ============================================================================

function interpolated_data = interpolate_data_band_specific_fixed(patient_data, band_params)
    % 🔧 JAVÍTOTT: Sáv-specifikus interpoláció proper datetime kezeléssel
    
    original_time = patient_data.Time;
    original_SE = patient_data.SE;
    original_rSO2 = patient_data.oper_side_oxig;
    original_MAP = patient_data.MAP;
    original_rSO2_contra = patient_data.other_side_oxig;
    
    % 🔥 DATETIME KONVERZIÓ (a működő alapkódból adaptálva)
    if isdatetime(original_time)
        % Datetime esetén másodpercekké konvertálás
        time_numeric = seconds(original_time - original_time(1));
    elseif isduration(original_time)
        % Duration esetén másodpercekké konvertálás
        time_numeric = seconds(original_time);
    else
        % Ha már szám, akkor használjuk úgy
        time_numeric = original_time;
    end
    
    % 🔥 ADAPTÍV IDŐABLAK KIVÁLASZTÁS
    total_duration = max(time_numeric) - min(time_numeric);
    target_duration = band_params.target_duration;
    
    if total_duration >= target_duration
        start_time = min(time_numeric) + (total_duration - target_duration) / 2;
        end_time = start_time + target_duration;
        time_mask = time_numeric >= start_time & time_numeric <= end_time;
        selected_time = time_numeric(time_mask);
        selected_SE = original_SE(time_mask);
        selected_rSO2 = original_rSO2(time_mask);
        selected_MAP = original_MAP(time_mask);
        selected_rSO2_contra = original_rSO2_contra(time_mask);
        actual_duration = target_duration;
    else
        selected_time = time_numeric;
        selected_SE = original_SE;
        selected_rSO2 = original_rSO2;
        selected_MAP = original_MAP;
        selected_rSO2_contra = original_rSO2_contra;
        actual_duration = total_duration;
    end
    
    % 🔥 VALID ADATOK SZŰRÉSE
    valid_mask = ~isnan(selected_time) & ~isnan(selected_SE) & ~isnan(selected_rSO2);
    
    if sum(valid_mask) < 3
        error('Elégtelen valid adatpontok az interpolációhoz (%d < 3)', sum(valid_mask));
    end
    
    valid_time = selected_time(valid_mask);
    valid_SE = selected_SE(valid_mask);
    valid_rSO2 = selected_rSO2(valid_mask);
    valid_MAP = selected_MAP(valid_mask);
    valid_rSO2_contra = selected_rSO2_contra(valid_mask);
    
    % 🔥 ÚJ IDŐRÁCS LÉTREHOZÁSA
    new_time = linspace(min(valid_time), max(valid_time), band_params.n_points);
    
    % 🔥 HELYES MINTAVÉTELI FREKVENCIA SZÁMÍTÁS
    actual_time_span = max(valid_time) - min(valid_time);
    dt = actual_time_span / (band_params.n_points - 1);
    actual_fs = 1 / dt;
    
    % 🔥 ROBUSZTUS INTERPOLÁCIÓ (mint a működő verzióban)
    interpolated_data = struct();
    try
        interpolated_data.SE = interp1(valid_time, valid_SE, new_time, 'pchip', 'extrap');
        interpolated_data.rSO2 = interp1(valid_time, valid_rSO2, new_time, 'pchip', 'extrap');
        interpolated_data.MAP = interp1(valid_time, valid_MAP, new_time, 'pchip', 'extrap');
        interpolated_data.rSO2_contra = interp1(valid_time, valid_rSO2_contra, new_time, 'pchip', 'extrap');
    catch
        % Fallback linear interpoláció
        interpolated_data.SE = interp1(valid_time, valid_SE, new_time, 'linear', 'extrap');
        interpolated_data.rSO2 = interp1(valid_time, valid_rSO2, new_time, 'linear', 'extrap');
        interpolated_data.MAP = interp1(valid_time, valid_MAP, new_time, 'linear', 'extrap');
        interpolated_data.rSO2_contra = interp1(valid_time, valid_rSO2_contra, new_time, 'linear', 'extrap');
    end
    
    % 🔥 JAVÍTOTT EREDMÉNYEK
    interpolated_data.time = new_time;
    interpolated_data.actual_duration = actual_time_span;  % Helyes időtartam
    interpolated_data.actual_fs = actual_fs;               % Helyes mintavételi frekvencia
    interpolated_data.dt = dt;                            % Helyes időlépés
    
    fprintf('    Interpoláció: %.1fs, %d pont, fs=%.3fHz, dt=%.3fs\n', ...
            actual_time_span, band_params.n_points, actual_fs, dt);
end

function [SE_filtered, rSO2_filtered] = apply_band_specific_filter_fixed(SE_signal, rSO2_signal, band_params, actual_fs)
    % 🔧 JAVÍTOTT: Sáv-specifikus bandpass szűrés HELYES actual_fs-sel
    
    if length(SE_signal) > 20 && actual_fs > 2 * band_params.range(2)  % Nyquist ellenőrzés
        try
            % Bandpass szűrő tervezése a HELYES mintavételi frekvenciával
            nyquist = actual_fs / 2;
            low_cutoff = max(band_params.range(1) / nyquist, 0.001);
            high_cutoff = min(band_params.range(2) / nyquist, 0.99);
            
            if low_cutoff < high_cutoff
                [b, a] = butter(4, [low_cutoff, high_cutoff], 'bandpass');
                SE_filtered = filtfilt(b, a, SE_signal);
                rSO2_filtered = filtfilt(b, a, rSO2_signal);
            else
                % Ha nem lehet bandpass, csak high-pass
                [b, a] = butter(2, low_cutoff, 'high');
                SE_filtered = filtfilt(b, a, SE_signal);
                rSO2_filtered = filtfilt(b, a, rSO2_signal);
            end
        catch
            % Ha szűrés nem sikerül, csak detrend
            SE_filtered = detrend(SE_signal);
            rSO2_filtered = detrend(rSO2_signal);
        end
    else
        SE_filtered = detrend(SE_signal);
        rSO2_filtered = detrend(rSO2_signal);
    end
end

%% KOMPREHENZÍV EXPORT FUNKCIÓ
% Ez a kód kiegészíti a főkódot - add hozzá a végére

%% 6. KOMPREHENZÍV BIOMARKER ÉS STATISZTIKAI EXPORT
fprintf('\n=== KOMPREHENZÍV EXPORT ===\n');

if n_successful_bands > 0
    % Comprehensive export minden sikeres sávra
    %export_comprehensive_biomarkers_and_stats(all_band_results, successful_bands);
    export_selected_angle_metrics(all_band_results, successful_bands);
    % Populációs összefoglaló export
    %export_population_summary(all_band_results, successful_bands);
else
    fprintf('❌ Nincs exportálható adat\n');
end

%% ============================================================================
%% KOMPREHENZÍV EXPORT FÜGGVÉNYEK
%% ============================================================================
function export_comprehensive_biomarkers_and_stats(all_band_results, successful_bands)
    % Komprehenzív biomarker és statisztikai export betegenként és sávonként
    
    fprintf('Komprehenzív biomarker export...\n');
    
    % 🔧 ELSŐ LÉPÉS: Összegyűjtjük az összes lehetséges oszlopnevet
    all_possible_columns = {};
    temp_rows = {};
    
    % Minden sáv minden betegének minden skálájának adatait összegyűjtjük
    row_counter = 1;
    for band_idx = 1:length(successful_bands)
        band_name = successful_bands{band_idx};
        band_data = all_band_results.(band_name);
        
        fprintf('  Feldolgozás: %s sáv (%d beteg)\n', band_name, length(band_data.patient_data));
        
        for patient_idx = 1:length(band_data.patient_data)
            patient_result = band_data.patient_data(patient_idx);
            
            % Minden skálára (scale_3, scale_5, scale_10, scale_20)
            scale_names = fieldnames(patient_result.SE_angle_analysis);
            
            for scale_idx = 1:length(scale_names)
                scale_name = scale_names{scale_idx};
                
                if isfield(patient_result.SE_angle_analysis, scale_name) && ...
                   isfield(patient_result.rSO2_angle_analysis, scale_name)
                    
                    % Egy sor adatot készítünk
                    row_data = extract_comprehensive_patient_data(patient_result, band_name, scale_name);
                    
                    % Tároljuk a sort és az oszlopneveket
                    temp_rows{row_counter} = row_data;
                    current_columns = row_data.Properties.VariableNames;
                    
                    % Újabb oszlopnevek hozzáadása
                    for col_idx = 1:length(current_columns)
                        col_name = current_columns{col_idx};
                        if ~ismember(col_name, all_possible_columns)
                            all_possible_columns{end+1} = col_name;
                        end
                    end
                    
                    row_counter = row_counter + 1;
                end
            end
        end
    end
    
    % 🔧 MÁSODIK LÉPÉS: Minden sort kiterjesztünk az összes oszloppal
    fprintf('  Oszlopok egységesítése (%d oszlop)...\n', length(all_possible_columns));
    
    standardized_rows = {};
    for i = 1:length(temp_rows)
        standardized_rows{i} = standardize_table_columns(temp_rows{i}, all_possible_columns);
    end
    
    % 🔧 HARMADIK LÉPÉS: Összes sor egyesítése
    if ~isempty(standardized_rows)
        all_patient_data = standardized_rows{1};
        for i = 2:length(standardized_rows)
            all_patient_data = [all_patient_data; standardized_rows{i}];
        end
        
        % CSV export
        export_to_csv(all_patient_data, 'comprehensive_biomarkers_and_stats.csv');
        fprintf('✅ Exportálva: comprehensive_biomarkers_and_stats.csv (%d sor, %d oszlop)\n', height(all_patient_data), width(all_patient_data));
    else
        fprintf('❌ Nincs exportálható adat\n');
    end
end

function standardized_table = standardize_table_columns(input_table, target_columns)
    % Táblázat oszlopainak egységesítése
    
    standardized_table = table();
    
    for i = 1:length(target_columns)
        col_name = target_columns{i};
        
        if ismember(col_name, input_table.Properties.VariableNames)
            % Oszlop létezik - átmásoljuk
            standardized_table.(col_name) = input_table.(col_name);
        else
            % Oszlop nem létezik - NaN-nal töltjük fel
            if height(input_table) > 0
                if contains(col_name, {'PatientID', 'SignalLength', 'N_Peaks'})
                    % Integer oszlopok
                    standardized_table.(col_name) = NaN(height(input_table), 1);
                elseif contains(col_name, {'Band', 'Scale'})
                    % String oszlopok
                    standardized_table.(col_name) = repmat({''}, height(input_table), 1);
                else
                    % Float oszlopok
                    standardized_table.(col_name) = NaN(height(input_table), 1);
                end
            end
        end
    end
end
function row_data = extract_comprehensive_patient_data(patient_result, band_name, scale_name)
    % Egyetlen beteg egyetlen sáv egyetlen skála adatainak kinyerése
    
    % Alap információk
    patient_id = patient_result.PatientID;
    actual_fs = patient_result.actual_fs;
    data_quality = patient_result.data_quality;
    signal_length = patient_result.signal_length;
    
    % SE és rSO2 angle analysis
    SE_analysis = patient_result.SE_angle_analysis.(scale_name);
    rSO2_analysis = patient_result.rSO2_angle_analysis.(scale_name);
     MAP_analysis = patient_result.MAP_angle_analysis.(scale_name);                      % ← ÚJ
    rSO2_contra_analysis = patient_result.rSO2_contra_angle_analysis.(scale_name);  
    
    % Cross-correlation analysis
    if isfield(patient_result.cross_angle_analysis, scale_name)
        cross_analysis = patient_result.cross_angle_analysis.(scale_name);
        cross_correlation = cross_analysis.correlation;
        circular_correlation = cross_analysis.circular_correlation;
        phase_lag = cross_analysis.phase_lag;
        max_cross_correlation = cross_analysis.max_cross_correlation;
    else
        cross_correlation = NaN;
        circular_correlation = NaN;
        phase_lag = NaN;
        max_cross_correlation = NaN;
    end
    
    % Band power
    SE_band_power = patient_result.SE_band_power;
    rSO2_band_power = patient_result.rSO2_band_power;
    MAP_band_power = patient_result.MAP_band_power;                      % ← ÚJ
    rSO2_contra_band_power = patient_result.rSO2_contra_band_power;  
    
    %% 🔥 RÉSZLETES SZÖGELOSZLÁS STATISZTIKÁK SZÁMÍTÁSA
    
    % SE szögeloszlás statisztikák
    SE_raw_angles = SE_analysis.raw_angles;
    [SE_q1, SE_q3] = calculate_quartiles(SE_raw_angles);
    SE_iqr = SE_q3 - SE_q1;
    SE_median = median(SE_raw_angles);
    SE_min = min(SE_raw_angles);
    SE_max = max(SE_raw_angles);
    SE_range = SE_max - SE_min;
    
    % rSO2 szögeloszlás statisztikák
    rSO2_raw_angles = rSO2_analysis.raw_angles;
    [rSO2_q1, rSO2_q3] = calculate_quartiles(rSO2_raw_angles);
    rSO2_iqr = rSO2_q3 - rSO2_q1;
    rSO2_median = median(rSO2_raw_angles);
    rSO2_min = min(rSO2_raw_angles);
    rSO2_max = max(rSO2_raw_angles);
    rSO2_range = rSO2_max - rSO2_min;

    MAP_raw_angles = MAP_analysis.raw_angles;
    [MAP_q1, MAP_q3] = calculate_quartiles(MAP_raw_angles);
    MAP_iqr = MAP_q3 - MAP_q1;
    MAP_median = median(MAP_raw_angles);
    MAP_min = min(MAP_raw_angles);
    MAP_max = max(MAP_raw_angles);
    MAP_range = MAP_max - MAP_min;
    
    % rSO2_contra szögeloszlás statisztikák                              % ← ÚJ BLOKK
    rSO2_contra_raw_angles = rSO2_contra_analysis.raw_angles;
    [rSO2_contra_q1, rSO2_contra_q3] = calculate_quartiles(rSO2_contra_raw_angles);
    rSO2_contra_iqr = rSO2_contra_q3 - rSO2_contra_q1;
    rSO2_contra_median = median(rSO2_contra_raw_angles);
    rSO2_contra_min = min(rSO2_contra_raw_angles);
    rSO2_contra_max = max(rSO2_contra_raw_angles);
    rSO2_contra_range = rSO2_contra_max - rSO2_contra_min;

    
    %% 🔧 BIZTONSÁGOS MEZŐ KIOLVASÁS
    % SE irány-specifikus statisztikák (biztonságos kiolvasás)
    SE_positive_mean = get_field_safe(SE_analysis, 'positive_mean', NaN);
    SE_positive_std = get_field_safe(SE_analysis, 'positive_std', NaN);
    SE_negative_mean = get_field_safe(SE_analysis, 'negative_mean', NaN);
    SE_negative_std = get_field_safe(SE_analysis, 'negative_std', NaN);
    
    % rSO2 irány-specifikus statisztikák (biztonságos kiolvasás)
    rSO2_positive_mean = get_field_safe(rSO2_analysis, 'positive_mean', NaN);
    rSO2_positive_std = get_field_safe(rSO2_analysis, 'positive_std', NaN);
    rSO2_negative_mean = get_field_safe(rSO2_analysis, 'negative_mean', NaN);
    rSO2_negative_std = get_field_safe(rSO2_analysis, 'negative_std', NaN);
    
    %% 🔥 BIOMARKEREK KINYERÉSE
    
    % SE biomarkerek
    if isfield(SE_analysis, 'band_biomarkers')
        SE_biomarkers = SE_analysis.band_biomarkers;
        SE_biomarker_names = fieldnames(SE_biomarkers);
        SE_biomarker_values = struct2array(SE_biomarkers);
    else
        SE_biomarker_names = {};
        SE_biomarker_values = [];
    end
    
    % rSO2 biomarkerek
    if isfield(rSO2_analysis, 'band_biomarkers')
        rSO2_biomarkers = rSO2_analysis.band_biomarkers;
        rSO2_biomarker_names = fieldnames(rSO2_biomarkers);
        rSO2_biomarker_values = struct2array(rSO2_biomarkers);
    else
        rSO2_biomarker_names = {};
        rSO2_biomarker_values = [];
    end
    
    %% 📊 TÁBLÁZAT SOR ÖSSZEÁLLÍTÁSA
    
    % Alap oszlopnevek és értékek
    row_data = table();
    
    % Egyszerű megközelítés: minden oszlop külön hozzáadása
    row_data.PatientID = patient_id;
    row_data.Band = {band_name};
    row_data.Scale = {scale_name};
    row_data.ActualFS = actual_fs;
    row_data.DataQuality = data_quality;
    row_data.SignalLength = signal_length;
    
    % SE alapstatisztikák
    row_data.SE_Mean = SE_analysis.mean_angle;
    row_data.SE_Std = SE_analysis.std_angle;
    row_data.SE_Median = SE_median;
    row_data.SE_Q1 = SE_q1;
    row_data.SE_Q3 = SE_q3;
    row_data.SE_IQR = SE_iqr;
    row_data.SE_Min = SE_min;
    row_data.SE_Max = SE_max;
    row_data.SE_Range = SE_range;
    row_data.SE_Kurtosis = SE_analysis.kurtosis;
    row_data.SE_Skewness = SE_analysis.skewness;
    
    % SE szögeloszlás metrikák
    row_data.SE_KDE_Entropy = SE_analysis.kde_entropy;
    row_data.SE_N_Peaks = SE_analysis.n_peaks;
    row_data.SE_Mode_Angle = SE_analysis.mode_angle;
    row_data.SE_Positive_Ratio = SE_analysis.positive_ratio;
    row_data.SE_Negative_Ratio = SE_analysis.negative_ratio;
    row_data.SE_Directional_Asymmetry = SE_analysis.directional_asymmetry;
    row_data.SE_Positive_Mean = SE_positive_mean;
    row_data.SE_Positive_Std = SE_positive_std;
    row_data.SE_Negative_Mean = SE_negative_mean;
    row_data.SE_Negative_Std = SE_negative_std;
    
    % rSO2 alapstatisztikák
    row_data.rSO2_Mean = rSO2_analysis.mean_angle;
    row_data.rSO2_Std = rSO2_analysis.std_angle;
    row_data.rSO2_Median = rSO2_median;
    row_data.rSO2_Q1 = rSO2_q1;
    row_data.rSO2_Q3 = rSO2_q3;
    row_data.rSO2_IQR = rSO2_iqr;
    row_data.rSO2_Min = rSO2_min;
    row_data.rSO2_Max = rSO2_max;
    row_data.rSO2_Range = rSO2_range;
    row_data.rSO2_Kurtosis = rSO2_analysis.kurtosis;
    row_data.rSO2_Skewness = rSO2_analysis.skewness;
    
    % rSO2 szögeloszlás metrikák
    row_data.rSO2_KDE_Entropy = rSO2_analysis.kde_entropy;
    row_data.rSO2_N_Peaks = rSO2_analysis.n_peaks;
    row_data.rSO2_Mode_Angle = rSO2_analysis.mode_angle;
    row_data.rSO2_Positive_Ratio = rSO2_analysis.positive_ratio;
    row_data.rSO2_Negative_Ratio = rSO2_analysis.negative_ratio;
    row_data.rSO2_Directional_Asymmetry = rSO2_analysis.directional_asymmetry;
    row_data.rSO2_Positive_Mean = rSO2_positive_mean;
    row_data.rSO2_Positive_Std = rSO2_positive_std;
    row_data.rSO2_Negative_Mean = rSO2_negative_mean;
    row_data.rSO2_Negative_Std = rSO2_negative_std;


    row_data.MAP_Mean = MAP_analysis.mean_angle;
    row_data.MAP_Std = MAP_analysis.std_angle;
    row_data.MAP_Median = MAP_median;
    row_data.MAP_Q1 = MAP_q1;
    row_data.MAP_Q3 = MAP_q3;
    row_data.MAP_IQR = MAP_iqr;
    row_data.MAP_Min = MAP_min;
    row_data.MAP_Max = MAP_max;
    row_data.MAP_Range = MAP_range;
    row_data.MAP_Kurtosis = MAP_analysis.kurtosis;
    row_data.MAP_Skewness = MAP_analysis.skewness;
    
    % MAP szögeloszlás metrikák                                          % ← ÚJ BLOKK
    row_data.MAP_KDE_Entropy = MAP_analysis.kde_entropy;
    row_data.MAP_N_Peaks = MAP_analysis.n_peaks;
    row_data.MAP_Mode_Angle = MAP_analysis.mode_angle;
    row_data.MAP_Positive_Ratio = MAP_analysis.positive_ratio;
    row_data.MAP_Negative_Ratio = MAP_analysis.negative_ratio;
    row_data.MAP_Directional_Asymmetry = MAP_analysis.directional_asymmetry;
    row_data.MAP_Positive_Mean = get_field_safe(MAP_analysis, 'positive_mean', NaN);
    row_data.MAP_Positive_Std = get_field_safe(MAP_analysis, 'positive_std', NaN);
    row_data.MAP_Negative_Mean = get_field_safe(MAP_analysis, 'negative_mean', NaN);
    row_data.MAP_Negative_Std = get_field_safe(MAP_analysis, 'negative_std', NaN);
    




    row_data.rSO2_contra_Mean = rSO2_contra_analysis.mean_angle;
    row_data.rSO2_contra_Std = rSO2_contra_analysis.std_angle;
    row_data.rSO2_contra_Median = rSO2_contra_median;
    row_data.rSO2_contra_Q1 = rSO2_contra_q1;
    row_data.rSO2_contra_Q3 = rSO2_contra_q3;
    row_data.rSO2_contra_IQR = rSO2_contra_iqr;
    row_data.rSO2_contra_Min = rSO2_contra_min;
    row_data.rSO2_contra_Max = rSO2_contra_max;
    row_data.rSO2_contra_Range = rSO2_contra_range;
    row_data.rSO2_contra_Kurtosis = rSO2_contra_analysis.kurtosis;
    row_data.rSO2_contra_Skewness = rSO2_contra_analysis.skewness;
    
    % rSO2_contra szögeloszlás metrikák                                          % ← ÚJ BLOKK
    row_data.rSO2_contra_KDE_Entropy = rSO2_contra_analysis.kde_entropy;
    row_data.rSO2_contra_N_Peaks = rSO2_contra_analysis.n_peaks;
    row_data.rSO2_contra_Mode_Angle = rSO2_contra_analysis.mode_angle;
    row_data.rSO2_contra_Positive_Ratio = rSO2_contra_analysis.positive_ratio;
    row_data.rSO2_contra_Negative_Ratio = rSO2_contra_analysis.negative_ratio;
    row_data.rSO2_contra_Directional_Asymmetry = rSO2_contra_analysis.directional_asymmetry;
    row_data.rSO2_contra_Positive_Mean = get_field_safe(rSO2_contra_analysis, 'positive_mean', NaN);
    row_data.rSO2_contra_Positive_Std = get_field_safe(rSO2_contra_analysis, 'positive_std', NaN);
    row_data.rSO2_contra_Negative_Mean = get_field_safe(rSO2_contra_analysis, 'negative_mean', NaN);
    row_data.rSO2_contra_Negative_Std = get_field_safe(rSO2_contra_analysis, 'negative_std', NaN);
    
    
    % Keresztkorreláció
    row_data.Cross_Correlation = cross_correlation;
    row_data.Circular_Correlation = circular_correlation;
    row_data.Phase_Lag = phase_lag;
    row_data.Max_Cross_Correlation = max_cross_correlation;
    
    % Spektrális teljesítmény
    row_data.SE_Band_Power = SE_band_power;
    row_data.rSO2_Band_Power = rSO2_band_power;
   
    row_data.MAP_Band_Power = MAP_band_power;                      % ← ÚJ
    row_data.rSO2_contra_Band_Power = rSO2_contra_band_power;
    
    % SE biomarkerek hozzáadása
    for i = 1:length(SE_biomarker_names)
        field_name = sprintf('SE_%s', SE_biomarker_names{i});
        row_data.(field_name) = SE_biomarker_values(i);
    end
    
    % rSO2 biomarkerek hozzáadása
    for i = 1:length(rSO2_biomarker_names)
        field_name = sprintf('rSO2_%s', rSO2_biomarker_names{i});
        row_data.(field_name) = rSO2_biomarker_values(i);
    end
end

function value = get_field_safe(struct_data, field_name, default_value)
    % Biztonságos mező kiolvasás
    if isfield(struct_data, field_name)
        value = struct_data.(field_name);
    else
        value = default_value;
    end
end
function [q1, q3] = calculate_quartiles(data)
    % Kvartilisek számítása
    sorted_data = sort(data);
    n = length(sorted_data);
    
    % Q1 (25. percentilis)
    q1_idx = 0.25 * (n + 1);
    if q1_idx == round(q1_idx)
        q1 = sorted_data(q1_idx);
    else
        lower_idx = floor(q1_idx);
        upper_idx = ceil(q1_idx);
        weight = q1_idx - lower_idx;
        q1 = sorted_data(lower_idx) * (1 - weight) + sorted_data(upper_idx) * weight;
    end
    
    % Q3 (75. percentilis)
    q3_idx = 0.75 * (n + 1);
    if q3_idx == round(q3_idx)
        q3 = sorted_data(q3_idx);
    else
        lower_idx = floor(q3_idx);
        upper_idx = ceil(q3_idx);
        weight = q3_idx - lower_idx;
        q3 = sorted_data(lower_idx) * (1 - weight) + sorted_data(upper_idx) * weight;
    end
end

function export_to_csv(data_table, filename)
    % Táblázat exportálása CSV-be
    try
        writetable(data_table, filename);
        fprintf('CSV exportálva: %s\n', filename);
    catch ME
        fprintf('❌ CSV export hiba: %s\n', ME.message);
    end
end

function export_population_summary(all_band_results, successful_bands)
    % Populációs szintű összefoglaló export
    
    fprintf('Populációs összefoglaló export...\n');
    
    summary_data = [];
    
    for band_idx = 1:length(successful_bands)
        band_name = successful_bands{band_idx};
        
        if isfield(all_band_results.(band_name), 'population_stats')
            pop_stats = all_band_results.(band_name).population_stats;
            
            % Egy sor a populációs statisztikákból
            row = table();
            row.Band = {band_name};
            row.N_Patients = pop_stats.n_patients;
            
            % SE populációs statisztikák
            row.SE_Mean_Angle_Pop_Mean = pop_stats.SE_mean_angle_mean;
            row.SE_Mean_Angle_Pop_Std = pop_stats.SE_mean_angle_std;
            row.SE_Entropy_Pop_Mean = pop_stats.SE_entropy_mean;
            row.SE_Entropy_Pop_Std = pop_stats.SE_entropy_std;
            row.SE_Asymmetry_Pop_Mean = pop_stats.SE_asymmetry_mean;
            row.SE_Asymmetry_Pop_Std = pop_stats.SE_asymmetry_std;
            
            % rSO2 populációs statisztikák
            row.rSO2_Mean_Angle_Pop_Mean = pop_stats.rSO2_mean_angle_mean;
            row.rSO2_Mean_Angle_Pop_Std = pop_stats.rSO2_mean_angle_std;
            row.rSO2_Entropy_Pop_Mean = pop_stats.rSO2_entropy_mean;
            row.rSO2_Entropy_Pop_Std = pop_stats.rSO2_entropy_std;
            row.rSO2_Asymmetry_Pop_Mean = pop_stats.rSO2_asymmetry_mean;
            row.rSO2_Asymmetry_Pop_Std = pop_stats.rSO2_asymmetry_std;
            
            summary_data = [summary_data; row];
        end
    end
    
    % Export
    if ~isempty(summary_data)
        export_to_csv(summary_data, 'population_summary_by_bands.csv');
        fprintf('✅ Populációs összefoglaló exportálva: population_summary_by_bands.csv\n');
    end
end

%% 📊 EXPORT PREVIEW FUNKCIÓ
function preview_export_structure(all_band_results, successful_bands)
    % Export struktúra előnézete
    
    fprintf('\n=== EXPORT STRUKTÚRA ELŐNÉZETE ===\n');
    
    if ~isempty(successful_bands)
        % Első sáv első betegének első skálája
        first_band = successful_bands{1};
        first_patient = all_band_results.(first_band).patient_data(1);
        first_scale = fieldnames(first_patient.SE_angle_analysis);
        first_scale = first_scale{1};
        
        % Példa sor generálása
        example_row = extract_comprehensive_patient_data(first_patient, first_band, first_scale);
        
        fprintf('Példa sor oszlopai (%d oszlop):\n', width(example_row));
        column_names = example_row.Properties.VariableNames;
        
        for i = 1:length(column_names)
            fprintf('%d. %s\n', i, column_names{i});
        end
        
        fprintf('\nPélda értékek (első 10 oszlop):\n');
        fprintf('PatientID: %d\n', example_row.PatientID);
        fprintf('Band: %s\n', example_row.Band{1});
        fprintf('Scale: %s\n', example_row.Scale{1});
        fprintf('SE_Mean: %.4f\n', example_row.SE_Mean);
        fprintf('SE_Std: %.4f\n', example_row.SE_Std);
        fprintf('SE_Q1: %.4f\n', example_row.SE_Q1);
        fprintf('SE_Q3: %.4f\n', example_row.SE_Q3);
        fprintf('SE_Kurtosis: %.4f\n', example_row.SE_Kurtosis);
        fprintf('SE_Skewness: %.4f\n', example_row.SE_Skewness);
        
        % Várható sorok száma becslése
        total_rows = 0;
        for i = 1:length(successful_bands)
            band_name = successful_bands{i};
            n_patients = length(all_band_results.(band_name).patient_data);
            n_scales = 4; % scale_3, scale_5, scale_10, scale_20
            total_rows = total_rows + n_patients * n_scales;
        end
        
        fprintf('\n📊 Várható export méret:\n');
        fprintf('Sávok: %d\n', length(successful_bands));
        fprintf('Összes sor: ~%d\n', total_rows);
        fprintf('Oszlopok: %d\n', width(example_row));
        fprintf('Becsült fájlméret: ~%.1f MB\n', total_rows * width(example_row) * 20 / 1024 / 1024);
    end
end

%% ============================================================================
%% KÖZÖS HELPER FÜGGVÉNYEK (A MŰKÖDŐ VERZIÓBÓL)
%% ============================================================================

function cleaned_signal = smart_gap_filling(signal, max_gap_points)
    cleaned_signal = signal;
    nan_positions = isnan(signal);
    
    if ~any(nan_positions)
        return;
    end
    
    gap_starts = find(diff([0; nan_positions]) == 1);
    gap_ends = find(diff([nan_positions; 0]) == -1);
    
    for i = 1:length(gap_starts)
        gap_length = gap_ends(i) - gap_starts(i) + 1;
        
        if gap_length <= max_gap_points
            start_idx = gap_starts(i);
            end_idx = gap_ends(i);
            
            before_val = NaN;
            after_val = NaN;
            
            if start_idx > 1
                before_val = signal(start_idx - 1);
            end
            if end_idx < length(signal)
                after_val = signal(end_idx + 1);
            end
            
            if ~isnan(before_val) && ~isnan(after_val)
                gap_values = linspace(before_val, after_val, gap_length + 2);
                cleaned_signal(start_idx:end_idx) = gap_values(2:end-1);
            elseif ~isnan(before_val)
                cleaned_signal(start_idx:end_idx) = before_val;
            elseif ~isnan(after_val)
                cleaned_signal(start_idx:end_idx) = after_val;
            end
        end
    end
end

function derivatives = calculate_multiscale_derivatives(signal, dt, window_sizes)
    % Multi-scale derivált számítás különböző módszerekkel
    derivatives = struct();
    
    % 1. Standard gradient
    derivatives.gradient = gradient(signal, dt);
    
    % 2. Különböző ablakméretű lokális deriváltak
    for i = 1:length(window_sizes)
        ws = window_sizes(i);
        scale_name = sprintf('scale_%d', ws);
        
        if ws <= length(signal)/4  % Csak ha van elég adat
            % Savitzky-Golay simított derivált
            if ws >= 3
                try
                    derivatives.(scale_name) = savitzky_golay_derivative(signal, ws, dt);
                catch
                    % Ha nincs Curve Fitting Toolbox, egyszerű finite difference
                    derivatives.(scale_name) = finite_difference_derivative(signal, ws, dt);
                end
            else
                derivatives.(scale_name) = gradient(signal, dt);
            end
        end
    end
end

function deriv = savitzky_golay_derivative(signal, window_size, dt)
    % Savitzky-Golay simított derivált
    if exist('sgolayfilt', 'file') == 2
        % Simítás előtt deriválás
        order = min(3, window_size-1);
        smoothed = sgolayfilt(signal, order, window_size);
        deriv = gradient(smoothed, dt);
    else
        % Fallback: egyszerű mozgóátlag + gradient
        smoothed = smoothdata(signal, 'movmean', window_size);
        deriv = gradient(smoothed, dt);
    end
end

function deriv = finite_difference_derivative(signal, window_size, dt)
    % Finite difference a megadott ablakmérettel
    deriv = zeros(size(signal));
    half_window = floor(window_size/2);
    
    for i = (half_window+1):(length(signal)-half_window)
        % Központi differencia az ablakon
        left_val = signal(i-half_window);
        right_val = signal(i+half_window);
        deriv(i) = (right_val - left_val) / (2 * half_window * dt);
    end
    
    % Széleken forward/backward difference
    deriv(1:half_window) = deriv(half_window+1);
    deriv((end-half_window+1):end) = deriv(end-half_window);
end

function [kde_values, kde_x] = calculate_kde(angles, bandwidth)
    % Kernel Density Estimation folytonos szögeloszláshoz
    kde_x = linspace(min(angles)-30, max(angles)+30, 500);
    kde_values = zeros(size(kde_x));
    
    % Gaussian kernel minden pontra
    for i = 1:length(angles)
        kernel_contrib = exp(-0.5 * ((kde_x - angles(i)) / bandwidth).^2);
        kde_values = kde_values + kernel_contrib;
    end
    
    % Normalizálás
    kde_values = kde_values / (length(angles) * bandwidth * sqrt(2*pi));
end

function analysis = analyze_continuous_angle_distribution(angles, kde_values, kde_x)
    % Folytonos szögeloszlás fejlett elemzése
    analysis = struct();
    
    % Alapstatisztikák
    analysis.mean_angle = mean(angles);
    analysis.std_angle = std(angles);
    analysis.kurtosis = kurtosis(angles);
    analysis.skewness = skewness(angles);
    analysis.range = max(angles) - min(angles);
    
    % KDE-alapú metrikák
    [~, max_idx] = max(kde_values);
    analysis.mode_angle = kde_x(max_idx);  % Módus
    
    % Többcsúcsúság detektálás
    [peaks, peak_locs] = findpeaks(kde_values, 'MinPeakHeight', max(kde_values)*0.1);
    analysis.n_peaks = length(peaks);
    analysis.peak_angles = kde_x(peak_locs);
    analysis.peak_heights = peaks;
    
    % Szóródási metrikák
    analysis.kde_entropy = -sum(kde_values(kde_values>0) .* log(kde_values(kde_values>0))) * (kde_x(2)-kde_x(1));
    
    % Irány-specifikus elemzés
    positive_angles = angles(angles > 0);  % Felfelé
    negative_angles = angles(angles < 0);  % Lefelé
    
    analysis.positive_ratio = length(positive_angles) / length(angles);
    analysis.negative_ratio = length(negative_angles) / length(angles);
    
    if ~isempty(positive_angles)
        analysis.positive_mean = mean(positive_angles);
        analysis.positive_std = std(positive_angles);
    end
    
    if ~isempty(negative_angles)
        analysis.negative_mean = mean(negative_angles);
        analysis.negative_std = std(negative_angles);
    end
    
    % Aszimmetria mérték
    analysis.directional_asymmetry = abs(analysis.positive_ratio - analysis.negative_ratio);
    
    % Tárolás a teljes eloszláshoz
    analysis.kde_x = kde_x;
    analysis.kde_values = kde_values;
    analysis.raw_angles = angles;
end

function wavelet_angles = calculate_wavelet_angle_distribution(signal, fs)
    % Wavelet-alapú frekvencia-függő szögeloszlás
    try
        % Continuous Wavelet Transform
        [wt, f] = cwt(signal, fs);
        
        % Instantaneous phase minden frekvencián
        phases = angle(wt);
        
        % Phase derivatives (frekvencia-függő irányek)
        phase_derivatives = diff(phases, 1, 2);  % Time derivative
        
        % Szögek konverziója
        wavelet_angles = struct();
        
        % Különböző frekvenciasávok
        freq_bands = [0.01 0.05; 0.05 0.1; 0.1 0.2; 0.2 0.4];
        band_names = {'VLF', 'LF', 'MF', 'HF'};
        
        for i = 1:size(freq_bands, 1)
            freq_mask = f >= freq_bands(i,1) & f <= freq_bands(i,2);
            if any(freq_mask)
                band_phases = phase_derivatives(freq_mask, :);
                band_angles = rad2deg(angle(mean(exp(1i*band_phases), 1)));
                wavelet_angles.(band_names{i}) = band_angles;
            end
        end
        
    catch
        wavelet_angles = [];
    end
end

function cross_analysis = analyze_cross_angle_correlations(SE_analysis, rSO2_analysis)
    % Keresztkorreláció SE és rSO2 szögeloszlások között
    cross_analysis = struct();
    
    scale_names = fieldnames(SE_analysis);
    
    for i = 1:length(scale_names)
        scale_name = scale_names{i};
        
        SE_angles = SE_analysis.(scale_name).raw_angles;
        rSO2_angles = rSO2_analysis.(scale_name).raw_angles;
        
        % Pearson korreláció
        [r, p] = corr(SE_angles', rSO2_angles');
        cross_analysis.(scale_name).correlation = r;
        cross_analysis.(scale_name).p_value = p;
        
        % Circular correlation (szögekhez megfelelőbb)
        circ_corr = circular_correlation(SE_angles, rSO2_angles);
        cross_analysis.(scale_name).circular_correlation = circ_corr;
        
        % Phase lag analysis
        [lag, max_corr] = analyze_phase_lag(SE_angles, rSO2_angles);
        cross_analysis.(scale_name).phase_lag = lag;
        cross_analysis.(scale_name).max_cross_correlation = max_corr;
    end
end

function circ_corr = circular_correlation(angles1, angles2)
    % Cirkuláris korreláció szögekhez
    % Konvertálás radiánba
    rad1 = deg2rad(angles1);
    rad2 = deg2rad(angles2);
    
    % Circular correlation coefficient
    sin1 = sin(rad1 - mean(rad1));
    sin2 = sin(rad2 - mean(rad2));
    
    circ_corr = sum(sin1 .* sin2) / sqrt(sum(sin1.^2) * sum(sin2.^2));
end

function [optimal_lag, max_correlation] = analyze_phase_lag(signal1, signal2)
    % Fázis késés elemzés keresztkorrelációval
    max_lag = min(50, floor(length(signal1)/4));  % Max lag
    
    [cross_corr, lags] = xcorr(signal1, signal2, max_lag, 'normalized');
    [max_correlation, max_idx] = max(abs(cross_corr));
    optimal_lag = lags(max_idx);
end

%% ============================================================================
%% SÁVÁSPECIFIKUS FÜGGVÉNYEK
%% ============================================================================

function biomarkers = calculate_band_biomarkers(angle_analysis, band_name)
    % Sáv-specifikus biomarkerek számítása
    biomarkers = struct();
    
    switch band_name
        case 'Endothelial'
            % Endothelial specifikus: metabolikus reguláció
            biomarkers.metabolic_regularity = 1 / (1 + angle_analysis.directional_asymmetry);
            biomarkers.endothelial_dysfunction_risk = (angle_analysis.directional_asymmetry > 0.3) * 1.0;
            biomarkers.no_mediated_capacity = min(1, angle_analysis.kde_entropy / 3.0);
            
        case 'Neurogenic'
            % Neurogenic specifikus: szimpatikus aktivitás
            biomarkers.sympathetic_activity = min(1, angle_analysis.kde_entropy / 4.0);
            biomarkers.baroreflex_sensitivity = 1 / (1 + abs(angle_analysis.mean_angle) / 10);
            biomarkers.autonomic_dysfunction_risk = (angle_analysis.kde_entropy < 2.5) * 1.0;
            
        case 'Myogenic'
            % Myogenic specifikus: autoregulációs kapacitás
            biomarkers.autoregulation_capacity = min(1, angle_analysis.n_peaks / 3);
            biomarkers.arterial_stiffness_index = angle_analysis.directional_asymmetry;
            biomarkers.bayliss_effect_integrity = (angle_analysis.n_peaks >= 2) * 1.0;
            
        otherwise
            biomarkers.general_regularity = 1 / (1 + angle_analysis.directional_asymmetry);
    end
end

function [SE_power, rSO2_power] = calculate_band_power(SE_signal, rSO2_signal, band_params, actual_fs)
    % Sáv-specifikus spektrális teljesítmény
    try
        % Power Spectral Density
        [SE_psd, f] = pwelch(SE_signal, [], [], [], actual_fs);
        [rSO2_psd, ~] = pwelch(rSO2_signal, [], [], [], actual_fs);
        
        % Sávban lévő teljesítmény
        band_mask = f >= band_params.range(1) & f <= band_params.range(2);
        
        if any(band_mask)
            SE_power = trapz(f(band_mask), SE_psd(band_mask));
            rSO2_power = trapz(f(band_mask), rSO2_psd(band_mask));
            
            % Relatív teljesítmény (teljes spektrumhoz képest)
            SE_power = SE_power / trapz(f, SE_psd);
            rSO2_power = rSO2_power / trapz(f, rSO2_psd);
        else
            SE_power = NaN;
            rSO2_power = NaN;
        end
        
    catch
        SE_power = NaN;
        rSO2_power = NaN;
    end
end

function pop_stats = calculate_band_population_statistics(patient_data_array, band_name)
    % Sáv-specifikus populációs statisztikák
    n_patients = length(patient_data_array);
    
    % Adatgyűjtés
    SE_mean_angles = zeros(n_patients, 1);
    rSO2_mean_angles = zeros(n_patients, 1);
    SE_entropies = zeros(n_patients, 1);
    rSO2_entropies = zeros(n_patients, 1);
    SE_asymmetries = zeros(n_patients, 1);
    rSO2_asymmetries = zeros(n_patients, 1);
    
    % Első skála adatainak használata (scale_3)
    for i = 1:n_patients
        patient = patient_data_array(i);
        
        if isfield(patient.SE_angle_analysis, 'scale_3')
            SE_analysis = patient.SE_angle_analysis.scale_3;
            rSO2_analysis = patient.rSO2_angle_analysis.scale_3;
            
            SE_mean_angles(i) = SE_analysis.mean_angle;
            rSO2_mean_angles(i) = rSO2_analysis.mean_angle;
            SE_entropies(i) = SE_analysis.kde_entropy;
            rSO2_entropies(i) = rSO2_analysis.kde_entropy;
            SE_asymmetries(i) = SE_analysis.directional_asymmetry;
            rSO2_asymmetries(i) = rSO2_analysis.directional_asymmetry;
        end
    end
    
    % Populációs statisztikák
    pop_stats = struct();
    pop_stats.band_name = band_name;
    pop_stats.n_patients = n_patients;
    
    pop_stats.SE_mean_angle_mean = mean(SE_mean_angles);
    pop_stats.SE_mean_angle_std = std(SE_mean_angles);
    pop_stats.rSO2_mean_angle_mean = mean(rSO2_mean_angles);
    pop_stats.rSO2_mean_angle_std = std(rSO2_mean_angles);
    
    pop_stats.SE_entropy_mean = mean(SE_entropies);
    pop_stats.SE_entropy_std = std(SE_entropies);
    pop_stats.rSO2_entropy_mean = mean(rSO2_entropies);
    pop_stats.rSO2_entropy_std = std(rSO2_entropies);
    
    pop_stats.SE_asymmetry_mean = mean(SE_asymmetries);
    pop_stats.SE_asymmetry_std = std(SE_asymmetries);
    pop_stats.rSO2_asymmetry_mean = mean(rSO2_asymmetries);
    pop_stats.rSO2_asymmetry_std = std(rSO2_asymmetries);
end

function visualize_band_specific_analysis(patient_result, window_sizes, band_name)
    % Sáv-specifikus vizualizáció
    patient_id = patient_result.PatientID;
    
    % 1. Eredeti jelek
    subplot(3,3,1);
    plot(patient_result.time_vector, patient_result.SE_signal, 'b-', 'LineWidth', 2);
    title(sprintf('%s SE - P%d', band_name, patient_id));
    xlabel('Idő (s)'); ylabel('SE');
    grid on;
    
    subplot(3,3,2);
    plot(patient_result.time_vector, patient_result.rSO2_signal, 'r-', 'LineWidth', 2);
    title(sprintf('%s rSO2 - P%d', band_name, patient_id));
    xlabel('Idő (s)'); ylabel('rSO2 (%)');
    grid on;
    
    % 2. Szögeloszlások az első skálán
    if isfield(patient_result.SE_angle_analysis, 'scale_3')
        subplot(3,3,3);
        SE_analysis = patient_result.SE_angle_analysis.scale_3;
        rSO2_analysis = patient_result.rSO2_angle_analysis.scale_3;
        
        plot(SE_analysis.kde_x, SE_analysis.kde_values, 'b-', 'LineWidth', 2, 'DisplayName', 'SE');
        hold on;
        plot(rSO2_analysis.kde_x, rSO2_analysis.kde_values, 'r-', 'LineWidth', 2, 'DisplayName', 'rSO2');
        xline(0, '--k');
        title(sprintf('%s Szögeloszlás', band_name));
        xlabel('Szög (°)'); ylabel('Sűrűség');
        legend('Location', 'best');
        grid on;
    end
    
    % 3. Keresztkorreláció
    subplot(3,3,4);
    if isfield(patient_result.cross_angle_analysis, 'scale_3')
        cross_corr = patient_result.cross_angle_analysis.scale_3.correlation;
        bar(1, cross_corr, 'FaceColor', [0.7 0.3 0.9]);
        title(sprintf('%s SE-rSO2 Korreláció', band_name));
        ylabel('r');
        ylim([-1 1]);
        grid on;
        set(gca, 'XTick', 1, 'XTickLabel', 'Scale-3');
    end
    
    % 4. Spektrális teljesítmény
    subplot(3,3,5);
    powers = [patient_result.SE_band_power, patient_result.rSO2_band_power];
    bar([1 2], powers, 'FaceColor', [0.3 0.7 0.5]);
    title(sprintf('%s Sáv Teljesítmény', band_name));
    set(gca, 'XTick', [1 2], 'XTickLabel', {'SE', 'rSO2'});
    ylabel('Relatív teljesítmény');
    grid on;
    
    % 5. Biomarkerek (ha vannak)
    if isfield(patient_result.SE_angle_analysis, 'scale_3') && ...
       isfield(patient_result.SE_angle_analysis.scale_3, 'band_biomarkers')
        subplot(3,3,6);
        biomarkers = patient_result.SE_angle_analysis.scale_3.band_biomarkers;
        marker_names = fieldnames(biomarkers);
        marker_values = struct2array(biomarkers);
        
        bar(1:length(marker_values), marker_values, 'FaceColor', [0.8 0.6 0.2]);
        title(sprintf('%s Biomarkerek', band_name));
        set(gca, 'XTick', 1:length(marker_names), 'XTickLabel', marker_names, 'XTickLabelRotation', 45);
        ylabel('Érték');
        grid on;
    end
    
    % 6. Információs panel
    subplot(3,3,7);
    text(0.1, 0.9, sprintf('Páciens: %d', patient_id), 'FontSize', 12);
    text(0.1, 0.8, sprintf('Sáv: %s', band_name), 'FontSize', 12);
    text(0.1, 0.7, sprintf('Fs: %.3f Hz', patient_result.actual_fs), 'FontSize', 12);
    text(0.1, 0.6, sprintf('Minőség: %.1f%%', patient_result.data_quality*100), 'FontSize', 12);
    text(0.1, 0.5, sprintf('Pontok: %d', patient_result.signal_length), 'FontSize', 12);
    axis off;
    title('Adatok');
    
    sgtitle(sprintf('%s Sáv Elemzés - Páciens %d (fs=%.3fHz)', band_name, patient_id, patient_result.actual_fs), ...
            'FontSize', 14, 'FontWeight', 'bold');
end

function visualize_multi_band_comparison(all_band_results, successful_bands)
    % Több sáv összehasonlító vizualizációja
    figure('Name', 'Multi-Sáv Összehasonlítás', 'Position', [100, 100, 1600, 1200]);
    
    n_bands = length(successful_bands);
    colors = lines(n_bands);
    
    % 1. Populációs átlagok összehasonlítása
    subplot(2,3,1);
    SE_means = zeros(n_bands, 1);
    SE_stds = zeros(n_bands, 1);
    rSO2_means = zeros(n_bands, 1);
    rSO2_stds = zeros(n_bands, 1);
    
    for i = 1:n_bands
        band_name = successful_bands{i};
        if isfield(all_band_results.(band_name), 'population_stats')
            stats = all_band_results.(band_name).population_stats;
            SE_means(i) = stats.SE_mean_angle_mean;
            SE_stds(i) = stats.SE_mean_angle_std;
            rSO2_means(i) = stats.rSO2_mean_angle_mean;
            rSO2_stds(i) = stats.rSO2_mean_angle_std;
        end
    end
    
    x = 1:n_bands;
    errorbar(x-0.15, SE_means, SE_stds, 'o-', 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'SE');
    hold on;
    errorbar(x+0.15, rSO2_means, rSO2_stds, 's-', 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'rSO2');
    
    set(gca, 'XTick', x, 'XTickLabel', successful_bands);
    title('Átlagos szögek sávonként');
    ylabel('Szög (°)');
    legend('Location', 'best');
    grid on;
    
    % 2. Sikerességi arányok
    subplot(2,3,2);
    success_rates = zeros(n_bands, 1);
    for i = 1:n_bands
        band_name = successful_bands{i};
        success_rates(i) = length(all_band_results.(band_name).patient_data);
    end
    
    bar(1:n_bands, success_rates, 'FaceColor', [0.6 0.8 0.4]);
    set(gca, 'XTick', 1:n_bands, 'XTickLabel', successful_bands);
    title('Sikeres páciensek száma');
    ylabel('Páciensek száma');
    grid on;
    
    % 3. Frekvenciatartományok
    subplot(2,3,3);
    for i = 1:n_bands
        band_name = successful_bands{i};
        band_params = all_band_results.(band_name).band_params;
        freq_range = band_params.range;
        
        plot([freq_range(1) freq_range(2)], [i i], 'o-', 'Color', colors(i,:), ...
             'LineWidth', 4, 'MarkerSize', 8, 'DisplayName', band_name);
        hold on;
    end
    
    set(gca, 'YTick', 1:n_bands, 'YTickLabel', successful_bands);
    xlabel('Frekvencia (Hz)');
    title('Frekvenciasávok');
    legend('Location', 'best');
    grid on;
    set(gca, 'XScale', 'log');
    
    % 4. Entrópia összehasonlítás
    subplot(2,3,4);
    for i = 1:n_bands
        band_name = successful_bands{i};
        if isfield(all_band_results.(band_name), 'population_stats')
            stats = all_band_results.(band_name).population_stats;
            SE_entropy = stats.SE_entropy_mean;
            rSO2_entropy = stats.rSO2_entropy_mean;
            
            bar(i-0.15, SE_entropy, 0.3, 'FaceColor', colors(1,:));
            hold on;
            bar(i+0.15, rSO2_entropy, 0.3, 'FaceColor', colors(2,:));
        end
    end
    
    set(gca, 'XTick', 1:n_bands, 'XTickLabel', successful_bands);
    title('KDE Entrópia sávonként');
    ylabel('Entrópia');
    legend({'SE', 'rSO2'}, 'Location', 'best');
    grid on;
    
    % 5. Aszimmetria összehasonlítás
    subplot(2,3,5);
    for i = 1:n_bands
        band_name = successful_bands{i};
        if isfield(all_band_results.(band_name), 'population_stats')
            stats = all_band_results.(band_name).population_stats;
            SE_asym = stats.SE_asymmetry_mean;
            rSO2_asym = stats.rSO2_asymmetry_mean;
            
            bar(i-0.15, SE_asym, 0.3, 'FaceColor', colors(1,:));
            hold on;
            bar(i+0.15, rSO2_asym, 0.3, 'FaceColor', colors(2,:));
        end
    end
    
    set(gca, 'XTick', 1:n_bands, 'XTickLabel', successful_bands);
    title('Irányaszimmetria sávonként');
    ylabel('Aszimmetria');
    legend({'SE', 'rSO2'}, 'Location', 'best');
    grid on;
    
    % 6. Sáv karakterisztikák
    subplot(2,3,6);
    text(0.1, 0.9, 'Sáv karakterisztikák:', 'FontSize', 14, 'FontWeight', 'bold');
    
    
    
    axis off;
    
    sgtitle('Multi-Sáv Összehasonlító Elemzés', 'FontSize', 16, 'FontWeight', 'bold');
end





function export_selected_angle_metrics(all_band_results, successful_bands)
    % Csak a kiválasztott metrikák exportálása
    % Formátum: sorok = metrikák, oszlopok = páciensek
    
    fprintf('\n=== KIVÁLASZTOTT METRIKÁK EXPORTÁLÁSA ===\n');
    
    % Kiválasztott metrikák definíciója
    selected_metrics = {
        'SE_IQR';
        'SE_Mode_Angle';
        'rSO2_Mean';
        'rSO2_Std';
        'rSO2_Max';
        'MAP_Negative_Ratio';
        'Phase_Lag';
        'MAP_Band_Power';
        'rSO2_baroreflex_sensitivity';
        'SE_bayliss_effect_integrity';
    };
    
    % Gyűjtsük össze az összes páciens ID-t
    all_patient_ids = [];
    for band_idx = 1:length(successful_bands)
        band_name = successful_bands{band_idx};
        band_data = all_band_results.(band_name);
        
        for patient_idx = 1:length(band_data.patient_data)
            patient_id = band_data.patient_data(patient_idx).PatientID;
            if ~ismember(patient_id, all_patient_ids)
                all_patient_ids = [all_patient_ids, patient_id];
            end
        end
    end
    
    n_patients = length(all_patient_ids);
    n_metrics = length(selected_metrics);
    
    % Adatmátrix inicializálása
    data_matrix = NaN(n_metrics, n_patients);
    
    % Minden pácienshez gyűjtsük ki a metrikákat
    for p_idx = 1:n_patients
        patient_id = all_patient_ids(p_idx);
        
        % Keressük meg ezt a pácienst a band results-ban
        for band_idx = 1:length(successful_bands)
            band_name = successful_bands{band_idx};
            band_data = all_band_results.(band_name);
            
            % Keressük meg a pácienst ebben a sávban
            patient_found = false;
            for patient_idx = 1:length(band_data.patient_data)
                if band_data.patient_data(patient_idx).PatientID == patient_id
                    patient_result = band_data.patient_data(patient_idx);
                    patient_found = true;
                    break;
                end
            end
            
            if patient_found
                % Használjuk az első elérhető skálát (általában scale_3)
                scale_names = fieldnames(patient_result.SE_angle_analysis);
                if ~isempty(scale_names)
                    scale_name = scale_names{1};  % Első skála
                    
                    % Metrikák kinyerése
                    metric_values = extract_selected_metrics(patient_result, band_name, scale_name);
                    
                    % Tároljuk az értékeket (felülírjuk NaN-okat ha van adat)
                    for m_idx = 1:n_metrics
                        if ~isnan(metric_values(m_idx))
                            data_matrix(m_idx, p_idx) = metric_values(m_idx);
                        end
                    end
                end
            end
        end
    end
    
    % CSV tábla létrehozása - verzió 1: sorok = metrikák
    csv_table = table(selected_metrics, 'VariableNames', {'Metric'});
    
    for p_idx = 1:n_patients
        patient_col_name = sprintf('Patient_%d', all_patient_ids(p_idx));
        csv_table.(patient_col_name) = data_matrix(:, p_idx);
    end
    
    % CSV mentése
    writetable(csv_table, 'selected_angle_metrics.csv');
    fprintf('✅ Exportálva: selected_angle_metrics.csv\n');
    fprintf('   Mátrix: %d metrika × %d páciens\n', n_metrics, n_patients);
    
    % Transzponált verzió - verzió 2: sorok = páciensek
    transposed_data = array2table(data_matrix', 'VariableNames', selected_metrics);
    transposed_data.PatientID = all_patient_ids';
    transposed_data = transposed_data(:, [end, 1:end-1]);  % PatientID első oszlopba
    
    writetable(transposed_data, 'selected_angle_metrics_transposed.csv');
    fprintf('✅ Exportálva: selected_angle_metrics_transposed.csv\n');
    fprintf('   Mátrix: %d páciens × %d metrika\n', n_patients, n_metrics);
    
    % Statisztikák
    fprintf('\n=== METRIKA STATISZTIKÁK ===\n');
    for m_idx = 1:n_metrics
        metric_values = data_matrix(m_idx, :);
        valid_values = metric_values(~isnan(metric_values));
        
        fprintf('\n%s:\n', selected_metrics{m_idx});
        fprintf('  Valid páciensek: %d/%d (%.1f%%)\n', ...
            length(valid_values), n_patients, length(valid_values)/n_patients*100);
        
        if ~isempty(valid_values)
            fprintf('  Átlag: %.4f ± %.4f\n', mean(valid_values), std(valid_values));
            fprintf('  Medián: %.4f\n', median(valid_values));
            fprintf('  Tartomány: [%.4f, %.4f]\n', min(valid_values), max(valid_values));
        end
    end
end

function metric_values = extract_selected_metrics(patient_result, band_name, scale_name)
    % Kiválasztott metrikák kinyerése egy pácienshez
    
    metric_values = NaN(10, 1);  % 10 metrika
    
    % Ellenőrizzük, hogy léteznek-e az elemzések
    if ~isfield(patient_result.SE_angle_analysis, scale_name) || ...
       ~isfield(patient_result.rSO2_angle_analysis, scale_name) || ...
       ~isfield(patient_result.MAP_angle_analysis, scale_name)
        return;
    end
    
    % Elemzések kinyerése
    SE_analysis = patient_result.SE_angle_analysis.(scale_name);
    rSO2_analysis = patient_result.rSO2_angle_analysis.(scale_name);
    MAP_analysis = patient_result.MAP_angle_analysis.(scale_name);
    
    % 1. SE_IQR
    if isfield(SE_analysis, 'raw_angles')
        SE_raw_angles = SE_analysis.raw_angles;
        SE_q1 = prctile(SE_raw_angles, 25);
        SE_q3 = prctile(SE_raw_angles, 75);
        metric_values(1) = SE_q3 - SE_q1;
    end
    
    % 2. SE_Mode_Angle
    if isfield(SE_analysis, 'mode_angle')
        metric_values(2) = SE_analysis.mode_angle;
    end
    
    % 3. rSO2_Mean
    if isfield(rSO2_analysis, 'mean_angle')
        metric_values(3) = rSO2_analysis.mean_angle;
    end
    
    % 4. rSO2_Std
    if isfield(rSO2_analysis, 'std_angle')
        metric_values(4) = rSO2_analysis.std_angle;
    end
    
    % 5. rSO2_Max
    if isfield(rSO2_analysis, 'raw_angles')
        metric_values(5) = max(rSO2_analysis.raw_angles);
    end
    
    % 6. MAP_Negative_Ratio
    if isfield(MAP_analysis, 'negative_ratio')
        metric_values(6) = MAP_analysis.negative_ratio;
    end
    
    % 7. Phase_Lag
    if isfield(patient_result.cross_angle_analysis, scale_name)
        cross_analysis = patient_result.cross_angle_analysis.(scale_name);
        if isfield(cross_analysis, 'phase_lag')
            metric_values(7) = cross_analysis.phase_lag;
        end
    end
    
    % 8. MAP_Band_Power
    if isfield(patient_result, 'MAP_band_power')
        metric_values(8) = patient_result.MAP_band_power;
    end
    
    % 9. rSO2_baroreflex_sensitivity (Neurogenic band specific)
    if strcmpi(band_name, 'Neurogenic')
        if isfield(rSO2_analysis, 'band_biomarkers') && ...
           isfield(rSO2_analysis.band_biomarkers, 'baroreflex_sensitivity')
            metric_values(9) = rSO2_analysis.band_biomarkers.baroreflex_sensitivity;
        end
    end
    
    % 10. SE_bayliss_effect_integrity (Myogenic band specific)
    if strcmpi(band_name, 'Myogenic')
        if isfield(SE_analysis, 'band_biomarkers') && ...
           isfield(SE_analysis.band_biomarkers, 'bayliss_effect_integrity')
            metric_values(10) = SE_analysis.band_biomarkers.bayliss_effect_integrity;
        end
    end
end
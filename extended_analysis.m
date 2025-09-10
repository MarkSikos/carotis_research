clear; clc;

%% AUTOREGULÁCIÓS ANALÍZIS - ADATFELKÉSZÍTÉS
% Csak interpoláció és preprocessing, metrikák nélkül

%% Paraméterek
data_file = 'data/df_unnormalized.csv';
target_oppart = 4;  % Clamp fázis
min_valid_ratio = 0.70;  % Min 70% valid adat kell

fprintf('=== AUTOREGULÁCIÓS ANALÍZIS - ADATFELKÉSZÍTÉS ===\n');

% Path hozzáadás
addpath(genpath('functions'));

%% Frekvencia sávok definíciója
freq_bands = struct();

% Endothelial - lassú, hosszú mérés
freq_bands.Endothelial = struct( ...
    'range', [0.003, 0.02], ...
    'optimal_duration', 1000, ...
    'target_fs', 0.1, ...
    'n_points', 100, ...
    'min_cycles', 3, ...
    'max_gap_seconds', 50, ...
    'detrend_method', 'emd', ...
    'artifact_z_threshold', 3.0, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

% Neurogenic - közepes
freq_bands.Neurogenic = struct( ...
    'range', [0.02, 0.06], ...
    'optimal_duration', 500, ...
    'target_fs', 0.2, ...
    'n_points', 100, ...
    'min_cycles', 4, ...
    'max_gap_seconds', 25, ...
    'detrend_method', 'emd', ...
    'artifact_z_threshold', 3.2, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

% Myogenic - gyors, rövid mérés
freq_bands.Myogenic = struct( ...
    'range', [0.06, 0.15], ...
    'optimal_duration', 300, ...
    'target_fs', 0.5, ...
    'n_points', 150, ...
    'min_cycles', 5, ...
    'max_gap_seconds', 12, ...
    'detrend_method', 'emd', ...
    'artifact_z_threshold', 3.5, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

% Respiratory - nagyon gyors
freq_bands.Respiratory = struct( ...
    'range', [0.15, 0.4], ...
    'optimal_duration', 200, ...
    'target_fs', 1.0, ...
    'n_points', 200, ...
    'min_cycles', 6, ...
    'max_gap_seconds', 6, ...
    'detrend_method', 'linear', ...
    'artifact_z_threshold', 4.0, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

% Cardiac - extrém gyors
freq_bands.Cardiac = struct( ...
    'range', [0.4, 2.0], ...
    'optimal_duration', 120, ...
    'target_fs', 5.0, ...
    'n_points', 600, ...
    'min_cycles', 10, ...
    'max_gap_seconds', 3, ...
    'detrend_method', 'linear', ...
    'artifact_z_threshold', 4.5, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

% VLF - Transfer Function
freq_bands.VLF = struct( ...
    'range', [0.02, 0.07], ...
    'optimal_duration', 600, ...
    'target_fs', 0.2, ...
    'n_points', 120, ...
    'min_cycles', 4, ...
    'max_gap_seconds', 30, ...
    'detrend_method', 'linear', ...
    'artifact_z_threshold', 3.0, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

% LF - Transfer Function
freq_bands.LF = struct( ...
    'range', [0.07, 0.2], ...
    'optimal_duration', 300, ...
    'target_fs', 0.5, ...
    'n_points', 150, ...
    'min_cycles', 5, ...
    'max_gap_seconds', 15, ...
    'detrend_method', 'linear', ...
    'artifact_z_threshold', 3.2, ...
    'physio_bounds', [0, 200; 0, 100; 30, 200]); % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max] - LAZÍTOTT!

band_names = fieldnames(freq_bands);

%% 1. Adatok betöltése
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



%% 2. SÁVONKÉNTI ADATFELKÉSZÍTÉS ÉS SPECIFIKUS INDEXEK
processed_data = struct();
specific_indices = struct(); % Új: specifikus indexek tárolása

for band_idx = 1:length(band_names)
    band_name = band_names{band_idx};
    band_params = freq_bands.(band_name);
    
    fprintf('\n=== %s SÁV ADATFELKÉSZÍTÉS ===\n', band_name);
    fprintf('Frekvencia: %.3f-%.3f Hz\n', band_params.range(1), band_params.range(2));
    fprintf('Target fs: %.3f Hz, Pontok: %d\n', band_params.target_fs, band_params.n_points);
    
    % Sáv-specifikus feldolgozott adatok tárolása
    band_data = struct();
    band_data.patient_data = [];
    band_data.valid_patients = [];
    
    %% Páciens loop - adatfelkészítés
    for p = 1:n_patients
        patient_id = patients(p);
        patient_data = phase_data(phase_data.Identifier == patient_id, :);
        
        % 1. Adatminőség ellenőrzés
        SE_signal = patient_data.SE;
        rSO2_signal = patient_data.oper_side_oxig;
        MAP_signal = patient_data.MAP;
        HR_signal = patient_data.HR;  % HR jel hozzáadása
        rSO2_contra_signal = patient_data.other_side_oxig;
        RES_signal = patient_data.resistance; % ÚJ: Resistance jel hozzáadása
        
        
        % Valid arány ellenőrzés (HR és Resistance is)
        se_valid_ratio = sum(~isnan(SE_signal)) / length(SE_signal);
        rso2_valid_ratio = sum(~isnan(rSO2_signal)) / length(rSO2_signal);
        map_valid_ratio = sum(~isnan(MAP_signal)) / length(MAP_signal);
        hr_valid_ratio = sum(~isnan(HR_signal)) / length(HR_signal);
        rso2_contra_valid_ratio = sum(~isnan(rSO2_contra_signal)) / length(rSO2_contra_signal);
        res_valid_ratio = sum(~isnan(RES_signal)) / length(RES_signal); % ÚJ: Resistance valid ratio

        if se_valid_ratio < min_valid_ratio || rso2_valid_ratio < min_valid_ratio || map_valid_ratio < min_valid_ratio || rso2_contra_valid_ratio < min_valid_ratio || res_valid_ratio < min_valid_ratio
            fprintf('  Páciens %d: elégtelen adatminőség (SE: %.1f%%, rSO2: %.1f%%, MAP: %.1f%%, HR: %.1f%%, rSO2_contra: %.1f%%, RES: %.1f%%)\n', ...
                    patient_id, se_valid_ratio*100, rso2_valid_ratio*100, map_valid_ratio*100, hr_valid_ratio*100, rso2_contra_valid_ratio*100, res_valid_ratio*100);
            continue;
        end
        
        % 2. Gap filling (HR és Resistance is)
        max_gap_points = round(band_params.max_gap_seconds * band_params.target_fs);

        SE_cleaned = smart_gap_filling(SE_signal, max_gap_points);
        rSO2_cleaned = smart_gap_filling(rSO2_signal, max_gap_points);
        MAP_cleaned = smart_gap_filling(MAP_signal, max_gap_points);
        HR_cleaned = smart_gap_filling(HR_signal, max_gap_points);
        rSO2_contra_cleaned = smart_gap_filling(rSO2_contra_signal, max_gap_points);
        RES_cleaned = smart_gap_filling(RES_signal, max_gap_points); % ÚJ: Resistance gap filling

        try
            % 3. Sáv-specifikus interpoláció
            patient_data_for_interp = patient_data;
            patient_data_for_interp.SE = SE_cleaned;
            patient_data_for_interp.oper_side_oxig = rSO2_cleaned;
            patient_data_for_interp.MAP = MAP_cleaned;
            patient_data_for_interp.HR = HR_cleaned;
            patient_data_for_interp.other_side_oxig = rSO2_contra_cleaned;
            patient_data_for_interp.resistance = RES_cleaned; % ÚJ: Resistance hozzáadása interpolációhoz
            interpolated_data = interpolate_data(patient_data_for_interp, ...
                                                band_params.n_points, ...
                                                band_params.optimal_duration);
            
            SE_interp = interpolated_data.SE;
            rSO2_interp = interpolated_data.rSO2;
            MAP_interp = interpolated_data.MAP;
            HR_interp = interpolated_data.HR;
            rSO2_contra_interp = interpolated_data.rSO2_contra;
            RES_interp = interpolated_data.resistance; % ÚJ: Resistance interpolált adat kinyerése
                        
            % 4. Missing values kezelése
            SE_clean = fillmissing(SE_interp, 'nearest');
            rSO2_clean = fillmissing(rSO2_interp, 'nearest');
            MAP_clean = fillmissing(MAP_interp, 'nearest');
            HR_clean = fillmissing(HR_interp, 'nearest');
            rSO2_contra_clean = fillmissing(rSO2_contra_interp, 'nearest');
            RES_clean = fillmissing(RES_interp, 'nearest'); % ÚJ: Resistance missing values kezelése
            
            % 4.5. FIZIOLÓGIAI HATÁROK ELLENŐRZÉSE (RAW ADATOKON!)
            physio_bounds = band_params.physio_bounds;
            SE_physio_valid = SE_clean >= physio_bounds(1,1) & SE_clean <= physio_bounds(1,2);
            rSO2_physio_valid = rSO2_clean >= physio_bounds(2,1) & rSO2_clean <= physio_bounds(2,2);
            MAP_physio_valid = MAP_clean >= physio_bounds(3,1) & MAP_clean <= physio_bounds(3,2);
            rSO2_contra_physio_valid = rSO2_contra_clean >= physio_bounds(2,1) & rSO2_contra_clean <= physio_bounds(2,2);
            %HR-nek nincs szigorú fiziológiai határ ellenőrzés (40-200 bpm tartomány túl széles)
            %Resistance-nek nincs szigorú fiziológiai határ ellenőrzés (értéktartomány változó lehet)
            
            % Valid arány RAW adatokon
            physio_valid_ratio = sum(SE_physio_valid & rSO2_physio_valid & MAP_physio_valid & rSO2_contra_physio_valid) / length(SE_clean);

            if physio_valid_ratio < min_valid_ratio
                fprintf('  Páciens %d: fiziológiai határok sértése (%.1f%% valid)\n', ...
                        patient_id, physio_valid_ratio*100);
                continue;
            end
            
            % 5. Detrending (sáv-specifikus) - CSAK ezután!
            SE_clean = advanced_detrend(SE_clean, band_params.detrend_method, 3);
            rSO2_clean = advanced_detrend(rSO2_clean, band_params.detrend_method, 3);
            MAP_clean = advanced_detrend(MAP_clean, band_params.detrend_method, 3);
            HR_clean = advanced_detrend(HR_clean, band_params.detrend_method, 3);
            rSO2_contra_clean = advanced_detrend(rSO2_contra_clean, band_params.detrend_method, 3);
            RES_clean = advanced_detrend(RES_clean, band_params.detrend_method, 3); % ÚJ: Resistance detrending

            % 6. Sáv-specifikus high-pass szűrés (egyszerűsített)
            cutoff_freq = max(0.001, band_params.range(1) * 0.5);
            if length(SE_clean) > 10 && cutoff_freq < band_params.target_fs/2
                % Egyszerű Butterworth high-pass szűrő (gyorsabb, stabilabb)
                try
                    [b, a] = butter(1, cutoff_freq/(band_params.target_fs/2), 'high');
                    SE_clean = filtfilt(b, a, SE_clean);
                    rSO2_clean = filtfilt(b, a, rSO2_clean);
                    MAP_clean = filtfilt(b, a, MAP_clean);
                    HR_clean = filtfilt(b, a, HR_clean);
                    rSO2_contra_clean = filtfilt(b, a, rSO2_contra_clean);
                    RES_clean = filtfilt(b, a, RES_clean); % ÚJ: Resistance szűrése
                catch
                    fprintf('    Szűrés kihagyva (nem kritikus)\n');
                end
            end
            
            % 7. Final validity check (feldolgozott adatok)
            % Most már csak finite értékeket ellenőrzünk
            SE_valid = isfinite(SE_clean);
            rSO2_valid = isfinite(rSO2_clean);
            MAP_valid = isfinite(MAP_clean);
            rSO2_contra_valid = isfinite(rSO2_contra_clean);
            RES_valid = isfinite(RES_clean); % ÚJ: Resistance validity check
            valid_ratio = sum(SE_valid & rSO2_valid & MAP_valid & rSO2_contra_valid & RES_valid) / length(SE_clean);

            if valid_ratio >= min_valid_ratio
                % Készítjük fel az adatokat a metrika számításokhoz
                processed_patient = struct();
                processed_patient.PatientID = patient_id;
                processed_patient.SE_clean = SE_clean;
                processed_patient.rSO2_clean = rSO2_clean;
                processed_patient.MAP_clean = MAP_clean;
                processed_patient.HR_clean = HR_clean;  % HR hozzáadása
                processed_patient.RES_clean = RES_clean; % ÚJ: Resistance hozzáadása
                processed_patient.valid_ratio = valid_ratio;
                processed_patient.sampling_rate = band_params.target_fs;
                processed_patient.freq_range = band_params.range;
                processed_patient.data_length_sec = length(SE_clean) / band_params.target_fs;
                
                % 🎯 SPECIFIKUS INDEXEK SZÁMÍTÁSA (sáv-függő)
                calculated_indices = calculate_band_specific_indices(SE_clean, rSO2_clean, HR_clean, rSO2_contra_clean, MAP_clean, RES_clean, band_name, band_params);
                % Indexek tárolása
                if ~isfield(specific_indices, band_name)
                    specific_indices.(band_name) = [];
                end
                
                % Páciens ID hozzáadása az indexekhez
                calculated_indices.PatientID = patient_id;
                specific_indices.(band_name) = [specific_indices.(band_name); calculated_indices];
                
                band_data.patient_data = [band_data.patient_data; processed_patient];
                band_data.valid_patients = [band_data.valid_patients; patient_id];
                
                fprintf('  Páciens %d: siker (%.1f%% valid, %.1f sec)\n', ...
                        patient_id, valid_ratio*100, processed_patient.data_length_sec);
            else
                fprintf('  Páciens %d: final validity check failed (%.1f%% valid)\n', ...
                        patient_id, valid_ratio*100);
            end
            
        catch ME
            fprintf('  Páciens %d: feldolgozási hiba - %s\n', patient_id, ME.message);
            continue;
        end
    end
    
    %% Sáv összegzés
    n_valid = length(band_data.valid_patients);
    fprintf('%s sáv: %d/%d páciens sikeresen feldolgozva (%.1f%%)\n', ...
            band_name, n_valid, n_patients, n_valid/n_patients*100);
    
    if n_valid >= 5
        % Statisztikák a feldolgozott adatokról
        all_SE = [];
        all_rSO2 = [];
        all_MAP = [];
        all_RES = []; % ÚJ: Resistance statisztikákhoz
        for i = 1:length(band_data.patient_data)
            all_SE = [all_SE; band_data.patient_data(i).SE_clean];
            all_rSO2 = [all_rSO2; band_data.patient_data(i).rSO2_clean];
            all_MAP = [all_MAP; band_data.patient_data(i).MAP_clean];
            all_RES = [all_RES; band_data.patient_data(i).RES_clean]; % ÚJ: Resistance adatok gyűjtése
        end
        
        fprintf('  SE átlag: %.1f ± %.1f mmHg\n', mean(all_SE), std(all_SE));
        fprintf('  rSO2 átlag: %.1f ± %.1f %%\n', mean(all_rSO2), std(all_rSO2));
        fprintf('  MAP átlag: %.1f ± %.1f mmHg\n', mean(all_MAP), std(all_MAP));
        fprintf('  Resistance átlag: %.1f ± %.1f\n', mean(all_RES), std(all_RES)); % ÚJ: Resistance statisztika kiírása
        fprintf('  Átlagos adathossz: %.1f sec\n', ...
                mean([band_data.patient_data.data_length_sec]));
        
        processed_data.(band_name) = band_data;
    else
        fprintf('  ELÉGTELEN ADAT - sáv kihagyva (minimum 5 páciens szükséges)\n');
    end
end

%% 3. Feldolgozott adatok mentése (opcionális)
fprintf('\n=== ADATFELKÉSZÍTÉS ÖSSZEGZÉS ===\n');
processed_bands = fieldnames(processed_data);
fprintf('Sikeresen feldolgozott sávok: %d/%d\n', length(processed_bands), length(band_names));

for i = 1:length(processed_bands)
    band_name = processed_bands{i};
    n_patients_band = length(processed_data.(band_name).valid_patients);
    fprintf('  %s: %d páciens\n', band_name, n_patients_band);
end

fprintf('\nAdatfelkészítés befejezve!\n');
fprintf('A processed_data változó tartalmazza az összes előkészített adatot.\n');

%% 4. SPECIFIKUS INDEXEK CSV EXPORTÁLÁSA
fprintf('\n=== SPECIFIKUS INDEXEK EXPORTÁLÁSA ===\n');
%export_specific_indices_to_csv(specific_indices, patients);

fprintf('\n=== SPECIFIKUS INDEXEK EXPORTÁLÁSA ===\n');
export_specific_indices_to_csv2(specific_indices, patients);

fprintf('Következő lépés: további metrika számítások hozzáadása.\n');
%% ========================================================================
%% SPECIFIKUS INDEXEK SZÁMÍTÁSA
%% ========================================================================
function indices = calculate_band_specific_indices(SE_clean, rSO2_clean, HR_clean, rSO2_contra_clean, MAP_clean, RES_clean, band_name, band_params)
% Sáv-specifikus indexek számítása - teljes metrika csomag
    
    indices = struct();
    fs = band_params.target_fs;
    
    switch lower(band_name)
        case 'endothelial'
            indices.autoregulation_effectiveness = calculate_autoregulation_effectiveness(SE_clean, rSO2_clean);
            indices.autoregulation_time_constant = calculate_autoregulation_time_constant(SE_clean, rSO2_clean, fs);
            bilateral_metrics = calculate_bilateral_autoregulation_metrics(rSO2_clean, rSO2_contra_clean, MAP_clean, fs);
            indices.bilateral_autoregulation_efficiency = bilateral_metrics.bilateral_autoregulation_efficiency;
            indices.HAI_COx = bilateral_metrics.HAI_COx;
            bootstrap_metrics = calculate_bootstrap_metrics(SE_clean, rSO2_clean);
            indices.COx_CI_lower = bootstrap_metrics.COx_CI(1);
            indices.COx_CI_upper = bootstrap_metrics.COx_CI(2);
            indices.COx_stability = bootstrap_metrics.stability;
            indices.cardiac_cerebral_coupling = calculate_cardiac_cerebral_coupling(HR_clean, rSO2_clean, MAP_clean, fs);
            % ÚJ: COx single hemisphere számítás (Endothelial sávban értelmes - 1000s duration)
            indices.COx_single_hemisphere = calculate_COx_single_hemisphere(rSO2_clean, MAP_clean, fs);
            % ÚJ: COx window számítás SE-rSO2 között (Endothelial sávban értelmes - 1000s duration)
            window_size = round(fs * 300); % 5 perces ablak
            indices.COx_window_SE_rSO2 = calculate_COx_window(SE_clean, rSO2_clean, window_size);
            % ÚJ: Cross-correlation comprehensive (Endothelial sávban értelmes - 1000s duration)
            xcorr_results = calculate_cross_correlation_comprehensive(SE_clean, rSO2_clean, fs);
            indices.xcorr_max_corr = xcorr_results.max_corr;
            indices.xcorr_optimal_lag = xcorr_results.optimal_lag;
            indices.xcorr_asymmetry = xcorr_results.asymmetry;
            indices.xcorr_width = xcorr_results.width;
            % ÚJ: Cross Sample Entropy (Endothelial sávban értelmes - 1000s duration)
            indices.cross_sampen = calculate_cross_sample_entropy(SE_clean, rSO2_clean);
            indices.CVRI = calculate_CVRI(RES_clean, MAP_clean, fs);
            dfa_results = calculate_DFA(SE_clean, rSO2_clean);
            indices.SE_alpha = dfa_results.SE_alpha;
            indices.rSO2_alpha = dfa_results.rSO2_alpha;
            indices.cross_DFA = dfa_results.cross_DFA;
            indices.estimated_CCP = calculate_estimated_CCP(SE_clean, rSO2_clean);
            indices.estimated_RAP = calculate_estimated_RAP(SE_clean, rSO2_clean);
            indices.fractal_dim_SE = calculate_fractal_dimension(SE_clean);
            indices.fractal_dim_rSO2 = calculate_fractal_dimension(rSO2_clean);
            hrx_results = calculate_heart_rate_reactivity_index(HR_clean, rSO2_clean, MAP_clean, fs);
            indices.HRx = hrx_results.HRx;
            indices.TOHRx = hrx_results.TOHRx;
            indices.hurst_SE = calculate_hurst_exponent(SE_clean);
            indices.hurst_rSO2 = calculate_hurst_exponent(rSO2_clean);
            mse_results = calculate_multiscale_entropy(SE_clean, rSO2_clean);
            indices.SE_MSE = mse_results.SE_MSE;
            indices.rSO2_MSE = mse_results.rSO2_MSE;
            indices.complexity_index = mse_results.complexity_index;
            indices.mutual_information = calculate_mutual_information(SE_clean, rSO2_clean);
            indices.Mx = calculate_Mx_comprehensive(SE_clean, rSO2_clean, fs);
            indices.normalized_mutual_information = calculate_normalized_mutual_information(SE_clean, rSO2_clean);
            indices.COx_temporal_variability = calculate_phase_COx_variability(rSO2_clean, MAP_clean, fs);
            indices.PRx = calculate_PRx_comprehensive(SE_clean, rSO2_clean, fs);
            indices.autoregulation_range = calculate_resistance_autoregulation_range(RES_clean, MAP_clean, fs);
            surrogate_results = calculate_surrogate_testing(SE_clean, rSO2_clean);
            indices.surrogate_p_value = surrogate_results.p_value;
            indices.surrogate_significance = surrogate_results.significance;

            endothelial_metrics = calculate_endothelial_specific_metrics(SE_clean, rSO2_clean, fs);
            indices.vasomotion_strength_SE = endothelial_metrics.vasomotion_strength_SE;
            indices.vasomotion_strength_rSO2 = endothelial_metrics.vasomotion_strength_rSO2;
            indices.vasomotion_coupling = endothelial_metrics.vasomotion_coupling;
            indices.endothelial_dysfunction_index = endothelial_metrics.endothelial_dysfunction_index;
            indices.vasomotion_regularity = endothelial_metrics.vasomotion_regularity;
            freq_results = calculate_frequency_domain_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.SE_power = freq_results.SE_power;
            indices.rSO2_power = freq_results.rSO2_power;
            indices.power_ratio = freq_results.power_ratio;
            indices.spectral_centroid_SE = freq_results.spectral_centroid_SE;
            indices.spectral_centroid_rSO2 = freq_results.spectral_centroid_rSO2;
            indices.bandwidth_SE = freq_results.bandwidth_SE;
            indices.bandwidth_rSO2 = freq_results.bandwidth_rSO2;
            indices.rolloff_SE = freq_results.rolloff_SE;
            indices.rolloff_rSO2 = freq_results.rolloff_rSO2;
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            indices.Lx = calculate_Lx_comprehensive(SE_clean, rSO2_clean, fs);
            % Endothelial case-hez hozzáadás (az autoregulation_range után):
            if std(RES_clean) > 0
                inverse_resistance = 1 ./ RES_clean;
                inverse_resistance(~isfinite(inverse_resistance)) = NaN;
                indices.FRx = calculate_COx_single_hemisphere(inverse_resistance, MAP_clean, fs);
            else
                indices.FRx = NaN;
            end
            if std(RES_clean) > 0 && std(rSO2_clean) > 0
                indices.vascular_compliance = calculate_vascular_compliance(RES_clean, rSO2_clean, MAP_clean, fs);
            else
                indices.vascular_compliance = NaN;
            end
            indices.Sx = calculate_Sx_comprehensive(SE_clean, rSO2_clean, fs);

            tf_results = calculate_time_frequency_analysis(SE_clean, rSO2_clean, fs, band_params.range);
            indices.phase_diff = tf_results.phase_diff;
            indices.amplitude_corr = tf_results.amplitude_corr;
            indices.phase_amp_coupling = tf_results.phase_amp_coupling;
            transfer_results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.transfer_gain = transfer_results.gain;
            indices.transfer_phase = transfer_results.phase;
            indices.transfer_coherence = transfer_results.coherence;
            indices.transfer_phase_lead = transfer_results.phase_lead;
            indices.transfer_gain_variability = transfer_results.gain_variability;
            indices.transfer_phase_variability = transfer_results.phase_variability;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            wavelet_results = calculate_wavelet_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.wavelet_coherence = wavelet_results.coherence;
            indices.wavelet_phase = wavelet_results.phase;
            indices.PLV = wavelet_results.PLV;
            indices.PPC = wavelet_results.PPC;
            indices.wPLI = wavelet_results.wPLI;
            fprintf('      %s: Wavelet Coherence = %.4f\n', band_name, indices.wavelet_coherence);
            fprintf('      %s: PLV = %.4f\n', band_name, indices.PLV);
            fprintf('      %s: wPLI = %.4f\n', band_name, indices.wPLI);
            fprintf('      %s: Transfer Gain = %.4f\n', band_name, indices.transfer_gain);
            fprintf('      %s: Transfer Phase = %.4f deg\n', band_name, indices.transfer_phase);
            fprintf('      %s: Transfer Coherence = %.4f\n', band_name, indices.transfer_coherence);
            fprintf('      %s: Phase Diff = %.4f\n', band_name, indices.phase_diff);
            fprintf('      %s: Amplitude Corr = %.4f\n', band_name, indices.amplitude_corr);
            fprintf('      %s: Phase-Amp Coupling = %.4f\n', band_name, indices.phase_amp_coupling);
            fprintf('      Endothelial: Sx = %.4f\n', indices.Sx);
            fprintf('      Endothelial: Vascular Compliance = %.4f\n', indices.vascular_compliance);
            fprintf('      Endothelial: FRx = %.4f\n', indices.FRx);
            fprintf('      %s: Lx = %.4f\n', band_name, indices.Lx);
            fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            fprintf('      Endothelial: SE Power = %.4f\n', indices.SE_power);
            fprintf('      Endothelial: Power Ratio = %.4f\n', indices.power_ratio);
            fprintf('      Endothelial: Vasomotion Strength SE = %.4f\n', indices.vasomotion_strength_SE);
            fprintf('      Endothelial: Vasomotion Strength rSO2 = %.4f\n', indices.vasomotion_strength_rSO2);
            fprintf('      Endothelial: Vasomotion Coupling = %.4f\n', indices.vasomotion_coupling);
            fprintf('      Endothelial: Endothelial Dysfunction Index = %.4f\n', indices.endothelial_dysfunction_index);
            fprintf('      Endothelial: Vasomotion Regularity = %.4f\n', indices.vasomotion_regularity);
            fprintf('      Endothelial: SE Alpha = %.4f\n', indices.SE_alpha);
            fprintf('      Endothelial: rSO2 Alpha = %.4f\n', indices.rSO2_alpha);
            fprintf('      Endothelial: Cross DFA = %.4f\n', indices.cross_DFA);
            fprintf('      Endothelial: CVRI = %.4f\n', indices.CVRI);
            fprintf('      Endothelial: Autoregulation Effectiveness = %.4f\n', indices.autoregulation_effectiveness);
            fprintf('      Endothelial: Autoregulation Time Constant = %.4f sec\n', indices.autoregulation_time_constant);
            fprintf('      Endothelial: Bilateral Efficiency = %.4f\n', indices.bilateral_autoregulation_efficiency);
            fprintf('      Endothelial: COx Stability = %.4f\n', indices.COx_stability);
            fprintf('      Endothelial: Cardiac-Cerebral Coupling = %.4f\n', indices.cardiac_cerebral_coupling);
            fprintf('      Endothelial: COx Single Hemisphere = %.4f\n', indices.COx_single_hemisphere);
            fprintf('      Endothelial: COx Window SE-rSO2 = %.4f\n', indices.COx_window_SE_rSO2);
            fprintf('      Endothelial: XCorr Max Corr = %.4f\n', indices.xcorr_max_corr);
            fprintf('      Endothelial: XCorr Optimal Lag = %.4f sec\n', indices.xcorr_optimal_lag);
            fprintf('      Endothelial: Cross SampEn = %.4f\n', indices.cross_sampen);
            fprintf('      Endothelial: Estimated CCP = %.4f\n', indices.estimated_CCP);
            fprintf('      Endothelial: Estimated RAP = %.4f\n', indices.estimated_RAP);
            fprintf('      Endothelial: Fractal Dim SE = %.4f\n', indices.fractal_dim_SE);
            fprintf('      Endothelial: Fractal Dim rSO2 = %.4f\n', indices.fractal_dim_rSO2);
            % cardiac_cerebral_coupling és autonomic_modulation_index már számítva van
            fprintf('      Endothelial: HRx = %.4f\n', indices.HRx);
            fprintf('      Endothelial: TOHRx = %.4f\n', indices.TOHRx);
            fprintf('      Endothelial: Hurst SE = %.4f\n', indices.hurst_SE);
            fprintf('      Endothelial: Hurst rSO2 = %.4f\n', indices.hurst_rSO2);
            fprintf('      Endothelial: SE MSE = %.4f\n', indices.SE_MSE);
            fprintf('      Endothelial: rSO2 MSE = %.4f\n', indices.rSO2_MSE);
            fprintf('      Endothelial: Complexity Index = %.4f\n', indices.complexity_index);
            fprintf('      Endothelial: Mutual Information = %.4f\n', indices.mutual_information);
            fprintf('      Endothelial: Mx = %.4f\n', indices.Mx);
            fprintf('      Endothelial: Normalized Mutual Information = %.4f\n', indices.normalized_mutual_information);
            fprintf('      Endothelial: COx Temporal Variability = %.4f\n', indices.COx_temporal_variability);
            fprintf('      Endothelial: PRx = %.4f\n', indices.PRx);
            fprintf('      Endothelial: Autoregulation Range = %.4f\n', indices.autoregulation_range);
            fprintf('      Endothelial: Surrogate P-value = %.4f\n', indices.surrogate_p_value);
            fprintf('      Endothelial: Surrogate Significance = %d\n', indices.surrogate_significance);


        case 'neurogenic'
            indices.autonomic_modulation_index = calculate_autonomic_modulation_index(HR_clean, rSO2_clean, fs);
            indices.autoregulation_effectiveness = calculate_autoregulation_effectiveness(SE_clean, rSO2_clean);
            indices.autoregulation_time_constant = calculate_autoregulation_time_constant(SE_clean, rSO2_clean, fs);
            bilateral_metrics = calculate_bilateral_autoregulation_metrics(rSO2_clean, rSO2_contra_clean, MAP_clean, fs);
            indices.bilateral_autoregulation_efficiency = bilateral_metrics.bilateral_autoregulation_efficiency;
            indices.HAI_COx = bilateral_metrics.HAI_COx;
            bootstrap_metrics = calculate_bootstrap_metrics(SE_clean, rSO2_clean);
            indices.COx_CI_lower = bootstrap_metrics.COx_CI(1);
            indices.COx_CI_upper = bootstrap_metrics.COx_CI(2);
            indices.COx_stability = bootstrap_metrics.stability;
            indices.cardiac_cerebral_coupling = calculate_cardiac_cerebral_coupling(HR_clean, rSO2_clean, MAP_clean, fs);
            % ÚJ: COx single hemisphere számítás (Neurogenic sávban értelmes - 500s duration)
            indices.COx_single_hemisphere = calculate_COx_single_hemisphere(rSO2_clean, MAP_clean, fs);
            % ÚJ: COx window számítás SE-rSO2 között (Neurogenic sávban értelmes - 500s duration)
            window_size = round(fs * 300); % 5 perces ablak
            indices.COx_window_SE_rSO2 = calculate_COx_window(SE_clean, rSO2_clean, window_size);
            % ÚJ: Cross-correlation comprehensive (Neurogenic sávban értelmes - 500s duration)
            xcorr_results = calculate_cross_correlation_comprehensive(SE_clean, rSO2_clean, fs);
            indices.xcorr_max_corr = xcorr_results.max_corr;
            indices.xcorr_optimal_lag = xcorr_results.optimal_lag;
            indices.xcorr_asymmetry = xcorr_results.asymmetry;
            indices.xcorr_width = xcorr_results.width;
            % ÚJ: Cross Sample Entropy (Neurogenic sávban értelmes - 500s duration)
            indices.cross_sampen = calculate_cross_sample_entropy(SE_clean, rSO2_clean);
            freq_results = calculate_frequency_domain_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.SE_power = freq_results.SE_power;
            indices.rSO2_power = freq_results.rSO2_power;
            indices.power_ratio = freq_results.power_ratio;
            indices.spectral_centroid_SE = freq_results.spectral_centroid_SE;
            indices.spectral_centroid_rSO2 = freq_results.spectral_centroid_rSO2;
            indices.bandwidth_SE = freq_results.bandwidth_SE;
            indices.bandwidth_rSO2 = freq_results.bandwidth_rSO2;
            indices.rolloff_SE = freq_results.rolloff_SE;
            indices.rolloff_rSO2 = freq_results.rolloff_rSO2;
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            indices.Lx = calculate_Lx_comprehensive(SE_clean, rSO2_clean, fs);
            neurogenic_metrics = calculate_neurogenic_specific_metrics(SE_clean, rSO2_clean, fs);
            indices.sympathetic_tone_SE = neurogenic_metrics.sympathetic_tone_SE;
            indices.sympathetic_tone_rSO2 = neurogenic_metrics.sympathetic_tone_rSO2;
            indices.autonomic_balance = neurogenic_metrics.autonomic_balance;
            indices.neurogenic_coupling_efficiency = neurogenic_metrics.neurogenic_coupling_efficiency;
            indices.neurogenic_burst_frequency = neurogenic_metrics.neurogenic_burst_frequency;
            indices.neurogenic_burst_duration = neurogenic_metrics.neurogenic_burst_duration;
            tf_results = calculate_time_frequency_analysis(SE_clean, rSO2_clean, fs, band_params.range);
            indices.phase_diff = tf_results.phase_diff;
            indices.amplitude_corr = tf_results.amplitude_corr;
            indices.phase_amp_coupling = tf_results.phase_amp_coupling;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            transfer_results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.transfer_gain = transfer_results.gain;
            indices.transfer_phase = transfer_results.phase;
            indices.transfer_coherence = transfer_results.coherence;
            indices.transfer_phase_lead = transfer_results.phase_lead;
            indices.transfer_gain_variability = transfer_results.gain_variability;
            indices.transfer_phase_variability = transfer_results.phase_variability;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            wavelet_results = calculate_wavelet_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.wavelet_coherence = wavelet_results.coherence;
            indices.wavelet_phase = wavelet_results.phase;
            indices.PLV = wavelet_results.PLV;
            indices.PPC = wavelet_results.PPC;
            indices.wPLI = wavelet_results.wPLI;
            fprintf('      %s: Wavelet Coherence = %.4f\n', band_name, indices.wavelet_coherence);
            fprintf('      %s: PLV = %.4f\n', band_name, indices.PLV);
            fprintf('      %s: wPLI = %.4f\n', band_name, indices.wPLI);
            fprintf('      %s: Transfer Gain = %.4f\n', band_name, indices.transfer_gain);
            fprintf('      %s: Transfer Phase = %.4f deg\n', band_name, indices.transfer_phase);
            fprintf('      %s: Transfer Coherence = %.4f\n', band_name, indices.transfer_coherence);
            fprintf('      %s: Phase Diff = %.4f\n', band_name, indices.phase_diff);
            fprintf('      %s: Amplitude Corr = %.4f\n', band_name, indices.amplitude_corr);
            fprintf('      %s: Phase-Amp Coupling = %.4f\n', band_name, indices.phase_amp_coupling);
            %fprintf('      Neurogenic: Sympathetic Tone SE = %.4f\n', indices.sympathetic_tone_SE);
            %fprintf('      Neurogenic: Sympathetic Tone rSO2 = %.4f\n', indices.sympathetic_tone_rSO2);
            %fprintf('      Neurogenic: Autonomic Balance = %.4f\n', indices.autonomic_balance);
            %fprintf('      Neurogenic: Neurogenic Coupling Efficiency = %.4f\n', indices.neurogenic_coupling_efficiency);
            %fprintf('      Neurogenic: Neurogenic Burst Frequency = %.4f\n', indices.neurogenic_burst_frequency);
            %fprintf('      Neurogenic: Neurogenic Burst Duration = %.4f\n', indices.neurogenic_burst_duration);
            %fprintf('      %s: Lx = %.4f\n', band_name, indices.Lx);
            %fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            %fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            %fprintf('      Neurogenic: SE Power = %.4f\n', indices.SE_power);
            %fprintf('      Neurogenic: Power Ratio = %.4f\n', indices.power_ratio);
            
            %fprintf('      Neurogenic: Autonomic Modulation Index = %.4f\n', indices.autonomic_modulation_index);
            %fprintf('      Neurogenic: Autoregulation Effectiveness = %.4f\n', indices.autoregulation_effectiveness);
            %fprintf('      Neurogenic: Autoregulation Time Constant = %.4f sec\n', indices.autoregulation_time_constant);
            %fprintf('      Neurogenic: Bilateral Efficiency = %.4f\n', indices.bilateral_autoregulation_efficiency);
            %fprintf('      Neurogenic: COx Stability = %.4f\n', indices.COx_stability);
            %fprintf('      Neurogenic: Cardiac-Cerebral Coupling = %.4f\n', indices.cardiac_cerebral_coupling);
            %fprintf('      Neurogenic: COx Single Hemisphere = %.4f\n', indices.COx_single_hemisphere);
            %fprintf('      Neurogenic: COx Window SE-rSO2 = %.4f\n', indices.COx_window_SE_rSO2);
            %fprintf('      Neurogenic: XCorr Max Corr = %.4f\n', indices.xcorr_max_corr);
            %fprintf('      Neurogenic: XCorr Optimal Lag = %.4f sec\n', indices.xcorr_optimal_lag);
            %fprintf('      Neurogenic: Cross SampEn = %.4f\n', indices.cross_sampen);
            
        case 'myogenic'
            bilateral_metrics = calculate_bilateral_autoregulation_metrics(rSO2_clean, rSO2_contra_clean, MAP_clean, fs);
            indices.bilateral_autoregulation_efficiency = bilateral_metrics.bilateral_autoregulation_efficiency;
            bootstrap_metrics = calculate_bootstrap_metrics(SE_clean, rSO2_clean);
            indices.COx_CI_lower = bootstrap_metrics.COx_CI(1);
            indices.COx_CI_upper = bootstrap_metrics.COx_CI(2);
            indices.COx_stability = bootstrap_metrics.stability;
            indices.cardiac_cerebral_coupling = calculate_cardiac_cerebral_coupling(HR_clean, rSO2_clean, MAP_clean, fs);
            % ÚJ: COx single hemisphere számítás (Myogenic sávban még éppen értelmes - 300s duration)
            indices.COx_single_hemisphere = calculate_COx_single_hemisphere(rSO2_clean, MAP_clean, fs);
            % ÚJ: COx window számítás SE-rSO2 között (Myogenic sávban még éppen értelmes - 300s duration)
            window_size = round(fs * 300); % 5 perces ablak
            indices.COx_window_SE_rSO2 = calculate_COx_window(SE_clean, rSO2_clean, window_size);
            % ÚJ: Cross-correlation comprehensive (Myogenic sávban határos - 300s duration)
            xcorr_results = calculate_cross_correlation_comprehensive(SE_clean, rSO2_clean, fs);
            indices.xcorr_max_corr = xcorr_results.max_corr;
            indices.xcorr_optimal_lag = xcorr_results.optimal_lag;
            indices.xcorr_asymmetry = xcorr_results.asymmetry;
            indices.xcorr_width = xcorr_results.width;
            % ÚJ: Cross Sample Entropy (Myogenic sávban határos - 300s duration)
            indices.cross_sampen = calculate_cross_sample_entropy(SE_clean, rSO2_clean);
            freq_results = calculate_frequency_domain_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.SE_power = freq_results.SE_power;
            indices.rSO2_power = freq_results.rSO2_power;
            indices.power_ratio = freq_results.power_ratio;
            indices.spectral_centroid_SE = freq_results.spectral_centroid_SE;
            indices.spectral_centroid_rSO2 = freq_results.spectral_centroid_rSO2;
            indices.bandwidth_SE = freq_results.bandwidth_SE;
            indices.bandwidth_rSO2 = freq_results.bandwidth_rSO2;
            indices.rolloff_SE = freq_results.rolloff_SE;
            indices.rolloff_rSO2 = freq_results.rolloff_rSO2;
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            indices.Lx = calculate_Lx_comprehensive(SE_clean, rSO2_clean, fs);
            myogenic_metrics = calculate_myogenic_specific_metrics(SE_clean, rSO2_clean, fs);
            indices.myogenic_reactivity = myogenic_metrics.myogenic_reactivity;
            indices.smooth_muscle_tone_SE = myogenic_metrics.smooth_muscle_tone_SE;
            indices.smooth_muscle_tone_rSO2 = myogenic_metrics.smooth_muscle_tone_rSO2;
            indices.myogenic_autoregulation_index = myogenic_metrics.myogenic_autoregulation_index;
            indices.vascular_compliance = myogenic_metrics.vascular_compliance;
            indices.myogenic_frequency_stability = myogenic_metrics.myogenic_frequency_stability;
            tf_results = calculate_time_frequency_analysis(SE_clean, rSO2_clean, fs, band_params.range);
            indices.phase_diff = tf_results.phase_diff;
            indices.amplitude_corr = tf_results.amplitude_corr;
            indices.phase_amp_coupling = tf_results.phase_amp_coupling;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            transfer_results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.transfer_gain = transfer_results.gain;
            indices.transfer_phase = transfer_results.phase;
            indices.transfer_coherence = transfer_results.coherence;
            indices.transfer_phase_lead = transfer_results.phase_lead;
            indices.transfer_gain_variability = transfer_results.gain_variability;
            indices.transfer_phase_variability = transfer_results.phase_variability;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            wavelet_results = calculate_wavelet_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.wavelet_coherence = wavelet_results.coherence;
            indices.wavelet_phase = wavelet_results.phase;
            indices.PLV = wavelet_results.PLV;
            indices.PPC = wavelet_results.PPC;
            indices.wPLI = wavelet_results.wPLI;
            fprintf('      %s: Wavelet Coherence = %.4f\n', band_name, indices.wavelet_coherence);
            fprintf('      %s: PLV = %.4f\n', band_name, indices.PLV);
            fprintf('      %s: wPLI = %.4f\n', band_name, indices.wPLI);
            fprintf('      %s: Transfer Gain = %.4f\n', band_name, indices.transfer_gain);
            fprintf('      %s: Transfer Phase = %.4f deg\n', band_name, indices.transfer_phase);
            fprintf('      %s: Transfer Coherence = %.4f\n', band_name, indices.transfer_coherence);
            fprintf('      %s: Phase Diff = %.4f\n', band_name, indices.phase_diff);
            fprintf('      %s: Amplitude Corr = %.4f\n', band_name, indices.amplitude_corr);
            fprintf('      %s: Phase-Amp Coupling = %.4f\n', band_name, indices.phase_amp_coupling);
            %fprintf('      Myogenic: Myogenic Reactivity = %.4f\n', indices.myogenic_reactivity);
            %fprintf('      Myogenic: Smooth Muscle Tone SE = %.4f\n', indices.smooth_muscle_tone_SE);
            %fprintf('      Myogenic: Smooth Muscle Tone rSO2 = %.4f\n', indices.smooth_muscle_tone_rSO2);
            %fprintf('      Myogenic: Myogenic Autoregulation Index = %.4f\n', indices.myogenic_autoregulation_index);
            %fprintf('      Myogenic: Vascular Compliance = %.4f\n', indices.vascular_compliance);
            %fprintf('      Myogenic: Myogenic Frequency Stability = %.4f\n', indices.myogenic_frequency_stability);
            %fprintf('      %s: Lx = %.4f\n', band_name, indices.Lx);
            %fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            %fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            %fprintf('      Myogenic: SE Power = %.4f\n', indices.SE_power);
            %fprintf('      Myogenic: Power Ratio = %.4f\n', indices.power_ratio);
            %fprintf('      Myogenic: Bilateral Efficiency = %.4f\n', indices.bilateral_autoregulation_efficiency);
            %fprintf('      Myogenic: COx Stability = %.4f\n', indices.COx_stability);
            %fprintf('      Myogenic: Cardiac-Cerebral Coupling = %.4f\n', indices.cardiac_cerebral_coupling);
            %fprintf('      Myogenic: COx Single Hemisphere = %.4f\n', indices.COx_single_hemisphere);
            %fprintf('      Myogenic: COx Window SE-rSO2 = %.4f\n', indices.COx_window_SE_rSO2);
            %fprintf('      Myogenic: XCorr Max Corr = %.4f\n', indices.xcorr_max_corr);
            %fprintf('      Myogenic: XCorr Optimal Lag = %.4f sec\n', indices.xcorr_optimal_lag);
            %fprintf('      Myogenic: Cross SampEn = %.4f\n', indices.cross_sampen);
            
        case 'respiratory'
            indices.autonomic_modulation_index = calculate_autonomic_modulation_index(HR_clean, rSO2_clean, fs);
            indices.cardiac_cerebral_coupling = calculate_cardiac_cerebral_coupling(HR_clean, rSO2_clean, MAP_clean, fs);
            % COx single hemisphere NINCS számítva (Respiratory sávban túl rövid - 200s duration)
            % XCorr NINCS számítva (Respiratory sávban túl rövid - 200s duration)
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            respiratory_metrics = calculate_respiratory_specific_metrics(SE_clean, rSO2_clean, fs);
            indices.respiratory_coupling_strength = respiratory_metrics.respiratory_coupling_strength;
            indices.respiratory_modulation_SE = respiratory_metrics.respiratory_modulation_SE;
            indices.respiratory_modulation_rSO2 = respiratory_metrics.respiratory_modulation_rSO2;
            indices.respiratory_phase_lag_sec = respiratory_metrics.respiratory_phase_lag_sec;
            indices.dominant_respiratory_frequency = respiratory_metrics.dominant_respiratory_frequency;
            indices.respiratory_rate_variability = respiratory_metrics.respiratory_rate_variability;
            fprintf('      Respiratory: Coupling Strength = %.4f\n', indices.respiratory_coupling_strength);
            fprintf('      Respiratory: Modulation SE = %.4f\n', indices.respiratory_modulation_SE);
            fprintf('      Respiratory: Modulation rSO2 = %.4f\n', indices.respiratory_modulation_rSO2);
            fprintf('      Respiratory: Phase Lag = %.4f sec\n', indices.respiratory_phase_lag_sec);
            fprintf('      Respiratory: Dominant Frequency = %.4f Hz\n', indices.dominant_respiratory_frequency);
            fprintf('      Respiratory: Rate Variability = %.4f\n', indices.respiratory_rate_variability);
            fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            fprintf('      Respiratory: Autonomic Modulation Index = %.4f\n', indices.autonomic_modulation_index);
            fprintf('      Respiratory: Cardiac-Cerebral Coupling = %.4f\n', indices.cardiac_cerebral_coupling);
            
        case 'cardiac'
            cardiac_metrics = calculate_cardiac_specific_metrics(SE_clean, rSO2_clean, fs);
            indices.cardiac_coupling_strength = cardiac_metrics.cardiac_coupling_strength;
            indices.pulse_amplitude_SE = cardiac_metrics.pulse_amplitude_SE;
            indices.pulse_amplitude_rSO2 = cardiac_metrics.pulse_amplitude_rSO2;
            indices.HRV_proxy = cardiac_metrics.HRV_proxy;
            indices.cardiac_phase_lag_sec = cardiac_metrics.cardiac_phase_lag_sec;
            indices.pulse_wave_velocity_proxy = cardiac_metrics.pulse_wave_velocity_proxy;
            indices.dominant_cardiac_frequency = cardiac_metrics.dominant_cardiac_frequency;
            indices.estimated_heart_rate = cardiac_metrics.estimated_heart_rate;
            % COx single hemisphere NINCS számítva (Cardiac sávban túl rövid - 120s duration)
            % XCorr NINCS számítva (Cardiac sávban túl rövid - 120s duration)
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            %fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            %fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            %fprintf('      Cardiac: Coupling Strength = %.4f\n', indices.cardiac_coupling_strength);
            %fprintf('      Cardiac: Pulse Amplitude SE = %.4f\n', indices.pulse_amplitude_SE);
            %fprintf('      Cardiac: HRV Proxy = %.4f\n', indices.HRV_proxy);
            %fprintf('      Cardiac: Estimated HR = %.1f bpm\n', indices.estimated_heart_rate);
            
        case 'vlf'
            indices.autoregulation_effectiveness = calculate_autoregulation_effectiveness(SE_clean, rSO2_clean);
            indices.autoregulation_time_constant = calculate_autoregulation_time_constant(SE_clean, rSO2_clean, fs);
            bilateral_metrics = calculate_bilateral_autoregulation_metrics(rSO2_clean, rSO2_contra_clean, MAP_clean, fs);
            indices.bilateral_autoregulation_efficiency = bilateral_metrics.bilateral_autoregulation_efficiency;
            indices.HAI_COx = bilateral_metrics.HAI_COx;
            bootstrap_metrics = calculate_bootstrap_metrics(SE_clean, rSO2_clean);
            indices.COx_CI_lower = bootstrap_metrics.COx_CI(1);
            indices.COx_CI_upper = bootstrap_metrics.COx_CI(2);
            indices.COx_stability = bootstrap_metrics.stability;
            indices.cardiac_cerebral_coupling = calculate_cardiac_cerebral_coupling(HR_clean, rSO2_clean, MAP_clean, fs);
            % ÚJ: COx single hemisphere számítás (VLF sávban értelmes - 600s duration)
            indices.COx_single_hemisphere = calculate_COx_single_hemisphere(rSO2_clean, MAP_clean, fs);
            % ÚJ: COx window számítás SE-rSO2 között (VLF sávban értelmes - 600s duration)
            window_size = round(fs * 300); % 5 perces ablak
            indices.COx_window_SE_rSO2 = calculate_COx_window(SE_clean, rSO2_clean, window_size);
            % ÚJ: Cross-correlation comprehensive (VLF sávban értelmes - 600s duration)
            xcorr_results = calculate_cross_correlation_comprehensive(SE_clean, rSO2_clean, fs);
            indices.xcorr_max_corr = xcorr_results.max_corr;
            indices.xcorr_optimal_lag = xcorr_results.optimal_lag;
            indices.xcorr_asymmetry = xcorr_results.asymmetry;
            indices.xcorr_width = xcorr_results.width;
            freq_results = calculate_frequency_domain_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.SE_power = freq_results.SE_power;
            indices.rSO2_power = freq_results.rSO2_power;
            indices.power_ratio = freq_results.power_ratio;
            indices.spectral_centroid_SE = freq_results.spectral_centroid_SE;
            indices.spectral_centroid_rSO2 = freq_results.spectral_centroid_rSO2;
            indices.bandwidth_SE = freq_results.bandwidth_SE;
            indices.bandwidth_rSO2 = freq_results.bandwidth_rSO2;
            indices.rolloff_SE = freq_results.rolloff_SE;
            indices.rolloff_rSO2 = freq_results.rolloff_rSO2;
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            indices.Lx = calculate_Lx_comprehensive(SE_clean, rSO2_clean, fs);
            tf_results = calculate_time_frequency_analysis(SE_clean, rSO2_clean, fs, band_params.range);
            indices.phase_diff = tf_results.phase_diff;
            indices.amplitude_corr = tf_results.amplitude_corr;
            indices.phase_amp_coupling = tf_results.phase_amp_coupling;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            transfer_results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.transfer_gain = transfer_results.gain;
            indices.transfer_phase = transfer_results.phase;
            indices.transfer_coherence = transfer_results.coherence;
            indices.transfer_phase_lead = transfer_results.phase_lead;
            indices.transfer_gain_variability = transfer_results.gain_variability;
            indices.transfer_phase_variability = transfer_results.phase_variability;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            wavelet_results = calculate_wavelet_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.wavelet_coherence = wavelet_results.coherence;
            indices.wavelet_phase = wavelet_results.phase;
            indices.PLV = wavelet_results.PLV;
            indices.PPC = wavelet_results.PPC;
            indices.wPLI = wavelet_results.wPLI;
            fprintf('      %s: Wavelet Coherence = %.4f\n', band_name, indices.wavelet_coherence);
            fprintf('      %s: PLV = %.4f\n', band_name, indices.PLV);
            fprintf('      %s: wPLI = %.4f\n', band_name, indices.wPLI);
            fprintf('      %s: Transfer Gain = %.4f\n', band_name, indices.transfer_gain);
            fprintf('      %s: Transfer Phase = %.4f deg\n', band_name, indices.transfer_phase);
            fprintf('      %s: Transfer Coherence = %.4f\n', band_name, indices.transfer_coherence);
            fprintf('      %s: Phase Diff = %.4f\n', band_name, indices.phase_diff);
            fprintf('      %s: Amplitude Corr = %.4f\n', band_name, indices.amplitude_corr);
            fprintf('      %s: Phase-Amp Coupling = %.4f\n', band_name, indices.phase_amp_coupling);
            % Calculate stability and classification metrics based on existing TF results
            if transfer_results.coherence >= 0.5
                if transfer_results.gain < 1.0
                    indices.autoregulation_status = 1; % Good autoregulation
                elseif transfer_results.gain >= 1.0 && transfer_results.gain < 1.5
                    indices.autoregulation_status = 0.5; % Impaired autoregulation
                else
                    indices.autoregulation_status = 0; % Poor autoregulation
                end
            else
                indices.autoregulation_status = NaN; % Unreliable due to low coherence
            end
            % Phase margin (stability measure)
            phase_margin_deg = 180 + transfer_results.phase;
            indices.phase_margin = phase_margin_deg;
            % Gain margin (stability measure)
            if transfer_results.gain > 0
                indices.gain_margin_db = -20 * log10(transfer_results.gain);
            else
                indices.gain_margin_db = NaN;
            end
            % Autoregulation efficiency
            indices.autoregulation_efficiency = transfer_results.coherence * (1 / (1 + transfer_results.gain));
            fprintf('      %s: Autoregulation Status = %.2f\n', band_name, indices.autoregulation_status);
            fprintf('      %s: Phase Margin = %.2f deg\n', band_name, indices.phase_margin);
            fprintf('      %s: Autoregulation Efficiency = %.4f\n', band_name, indices.autoregulation_efficiency);
            %fprintf('      %s: Lx = %.4f\n', band_name, indices.Lx);
            %fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            %fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            %fprintf('      VLF: SE Power = %.4f\n', indices.SE_power);
            %fprintf('      VLF: Power Ratio = %.4f\n', indices.power_ratio);
            %fprintf('      VLF: Autoregulation Effectiveness = %.4f\n', indices.autoregulation_effectiveness);
            %fprintf('      VLF: Autoregulation Time Constant = %.4f sec\n', indices.autoregulation_time_constant);
            %fprintf('      VLF: Bilateral Efficiency = %.4f\n', indices.bilateral_autoregulation_efficiency);
            %fprintf('      VLF: COx Stability = %.4f\n', indices.COx_stability);
            %fprintf('      VLF: Cardiac-Cerebral Coupling = %.4f\n', indices.cardiac_cerebral_coupling);
            %fprintf('      VLF: COx Single Hemisphere = %.4f\n', indices.COx_single_hemisphere);
            %fprintf('      VLF: COx Window SE-rSO2 = %.4f\n', indices.COx_window_SE_rSO2);
            %fprintf('      VLF: XCorr Max Corr = %.4f\n', indices.xcorr_max_corr);
            %fprintf('      VLF: XCorr Optimal Lag = %.4f sec\n', indices.xcorr_optimal_lag);
            
        case 'lf'
            % ÚJ: COx single hemisphere számítás (LF sávban még éppen értelmes - 300s duration)
            indices.COx_single_hemisphere = calculate_COx_single_hemisphere(rSO2_clean, MAP_clean, fs);
            % ÚJ: COx window számítás SE-rSO2 között (LF sávban még éppen értelmes - 300s duration)
            window_size = round(fs * 300); % 5 perces ablak
            indices.COx_window_SE_rSO2 = calculate_COx_window(SE_clean, rSO2_clean, window_size);
            % ÚJ: Cross-correlation comprehensive (LF sávban határos - 300s duration)
            xcorr_results = calculate_cross_correlation_comprehensive(SE_clean, rSO2_clean, fs);
            indices.xcorr_max_corr = xcorr_results.max_corr;
            indices.xcorr_optimal_lag = xcorr_results.optimal_lag;
            indices.xcorr_asymmetry = xcorr_results.asymmetry;
            indices.xcorr_width = xcorr_results.width;
            freq_results = calculate_frequency_domain_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.SE_power = freq_results.SE_power;
            indices.rSO2_power = freq_results.rSO2_power;
            indices.power_ratio = freq_results.power_ratio;
            indices.spectral_centroid_SE = freq_results.spectral_centroid_SE;
            indices.spectral_centroid_rSO2 = freq_results.spectral_centroid_rSO2;
            indices.bandwidth_SE = freq_results.bandwidth_SE;
            indices.bandwidth_rSO2 = freq_results.bandwidth_rSO2;
            indices.rolloff_SE = freq_results.rolloff_SE;
            indices.rolloff_rSO2 = freq_results.rolloff_rSO2;
            granger_results = calculate_granger_causality(SE_clean, rSO2_clean);
            indices.SE_to_rSO2_causality = granger_results.SE_to_rSO2;
            indices.rSO2_to_SE_causality = granger_results.rSO2_to_SE;
            indices.bidirectional_causality = granger_results.bidirectional;
            indices.Lx = calculate_Lx_comprehensive(SE_clean, rSO2_clean, fs);
            tf_results = calculate_time_frequency_analysis(SE_clean, rSO2_clean, fs, band_params.range);
            indices.phase_diff = tf_results.phase_diff;
            indices.amplitude_corr = tf_results.amplitude_corr;
            indices.phase_amp_coupling = tf_results.phase_amp_coupling;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            transfer_results = calculate_transfer_function_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.transfer_gain = transfer_results.gain;
            indices.transfer_phase = transfer_results.phase;
            indices.transfer_coherence = transfer_results.coherence;
            indices.transfer_phase_lead = transfer_results.phase_lead;
            indices.transfer_gain_variability = transfer_results.gain_variability;
            indices.transfer_phase_variability = transfer_results.phase_variability;
            % Endothelial, Neurogenic, Myogenic, VLF, LF case-ekbe:
            wavelet_results = calculate_wavelet_comprehensive(SE_clean, rSO2_clean, fs, band_params.range);
            indices.wavelet_coherence = wavelet_results.coherence;
            indices.wavelet_phase = wavelet_results.phase;
            indices.PLV = wavelet_results.PLV;
            indices.PPC = wavelet_results.PPC;
            indices.wPLI = wavelet_results.wPLI;
            fprintf('      %s: Wavelet Coherence = %.4f\n', band_name, indices.wavelet_coherence);
            fprintf('      %s: PLV = %.4f\n', band_name, indices.PLV);
            fprintf('      %s: wPLI = %.4f\n', band_name, indices.wPLI);
            fprintf('      %s: Transfer Gain = %.4f\n', band_name, indices.transfer_gain);
            fprintf('      %s: Transfer Phase = %.4f deg\n', band_name, indices.transfer_phase);
            fprintf('      %s: Transfer Coherence = %.4f\n', band_name, indices.transfer_coherence);
            fprintf('      %s: Phase Diff = %.4f\n', band_name, indices.phase_diff);
            fprintf('      %s: Amplitude Corr = %.4f\n', band_name, indices.amplitude_corr);
            fprintf('      %s: Phase-Amp Coupling = %.4f\n', band_name, indices.phase_amp_coupling);
            % Calculate stability and classification metrics based on existing TF results
            if transfer_results.coherence >= 0.5
                if transfer_results.gain < 1.0
                    indices.autoregulation_status = 1; % Good autoregulation
                elseif transfer_results.gain >= 1.0 && transfer_results.gain < 1.5
                    indices.autoregulation_status = 0.5; % Impaired autoregulation
                else
                    indices.autoregulation_status = 0; % Poor autoregulation
                end
            else
                indices.autoregulation_status = NaN; % Unreliable due to low coherence
            end
            % Phase margin (stability measure)
            phase_margin_deg = 180 + transfer_results.phase;
            indices.phase_margin = phase_margin_deg;
            % Gain margin (stability measure)
            if transfer_results.gain > 0
                indices.gain_margin_db = -20 * log10(transfer_results.gain);
            else
                indices.gain_margin_db = NaN;
            end
            % Autoregulation efficiency
            indices.autoregulation_efficiency = transfer_results.coherence * (1 / (1 + transfer_results.gain));
            fprintf('      %s: Autoregulation Status = %.2f\n', band_name, indices.autoregulation_status);
            fprintf('      %s: Phase Margin = %.2f deg\n', band_name, indices.phase_margin);
            fprintf('      %s: Autoregulation Efficiency = %.4f\n', band_name, indices.autoregulation_efficiency);
            
            %fprintf('      %s: Lx = %.4f\n', band_name, indices.Lx);
            %fprintf('      %s: SE->rSO2 Causality = %.4f\n', band_name, indices.SE_to_rSO2_causality);
            %fprintf('      %s: rSO2->SE Causality = %.4f\n', band_name, indices.rSO2_to_SE_causality);
            %fprintf('      LF: SE Power = %.4f\n', indices.SE_power);
            %fprintf('      LF: Power Ratio = %.4f\n', indices.power_ratio);
            %fprintf('      LF: COx Single Hemisphere = %.4f\n', indices.COx_single_hemisphere);
            %fprintf('      LF: COx Window SE-rSO2 = %.4f\n', indices.COx_window_SE_rSO2);
            %fprintf('      LF: XCorr Max Corr = %.4f\n', indices.xcorr_max_corr);
            %fprintf('      LF: XCorr Optimal Lag = %.4f sec\n', indices.xcorr_optimal_lag);
            
        otherwise
            fprintf('      %s: Unknown band\n', band_name);
    end
end

function export_specific_indices_to_csv(specific_indices, all_patients)
    % CSV export: X dim = indexek, Y dim = páciensek
    
    % Sáv-index mátrix definíciója - CARDIAC-SPECIFIC METRIKÁK + COx metrikák + ÚJ XCorr metrikák HOZZÁADVA
    band_index_map = containers.Map();
    band_index_map('Endothelial') = {'autoregulation_effectiveness', 'autoregulation_time_constant', 'bilateral_autoregulation_efficiency', 'HAI_COx', 'COx_CI_lower', 'COx_CI_upper', 'COx_stability', 'cardiac_cerebral_coupling', 'COx_single_hemisphere', 'COx_window_SE_rSO2', 'xcorr_max_corr', 'xcorr_optimal_lag', 'xcorr_asymmetry', 'xcorr_width', 'CVRI', 'SE_alpha', 'rSO2_alpha', 'cross_DFA','vasomotion_strength_SE', 'vasomotion_strength_rSO2', 'vasomotion_coupling', 'endothelial_dysfunction_index', 'vasomotion_regularity','estimated_CCP', 'estimated_RAP','fractal_dim_SE', 'fractal_dim_rSO2','SE_power', 'rSO2_power', 'power_ratio', 'spectral_centroid_SE', 'spectral_centroid_rSO2', 'bandwidth_SE', 'bandwidth_rSO2', 'rolloff_SE', 'rolloff_rSO2','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality','HRx', 'TOHRx','hurst_SE', 'hurst_rSO2','Lx','SE_MSE', 'rSO2_MSE', 'complexity_index','mutual_information','Mx', 'normalized_mutual_information','COx_temporal_variability','PRx', 'autoregulation_range','FRx','vascular_compliance','surrogate_p_value','surrogate_significance','Sx','phase_diff','amplitude_corr','phase_amp_coupling', 'transfer_gain','transfer_phase','transfer_coherence','transfer_phase_lead','transfer_gain_variability','transfer_phase_variability','wavelet_coherence','wavelet_phase','PLV','PPC','wPLI'};
    band_index_map('Neurogenic') = {'autonomic_modulation_index', 'autoregulation_effectiveness', 'autoregulation_time_constant', 'bilateral_autoregulation_efficiency', 'HAI_COx', 'COx_CI_lower', 'COx_CI_upper', 'COx_stability', 'cardiac_cerebral_coupling', 'COx_single_hemisphere', 'COx_window_SE_rSO2', 'xcorr_max_corr', 'xcorr_optimal_lag', 'xcorr_asymmetry', 'xcorr_width','SE_power', 'rSO2_power', 'power_ratio', 'spectral_centroid_SE', 'spectral_centroid_rSO2', 'bandwidth_SE', 'bandwidth_rSO2', 'rolloff_SE', 'rolloff_rSO2','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality','Lx','sympathetic_tone_SE', 'sympathetic_tone_rSO2', 'autonomic_balance', 'neurogenic_coupling_efficiency', 'neurogenic_burst_frequency', 'neurogenic_burst_duration','phase_diff','amplitude_corr','phase_amp_coupling', 'transfer_gain','transfer_phase','transfer_coherence','transfer_phase_lead','transfer_gain_variability','transfer_phase_variability','wavelet_coherence','wavelet_phase','PLV','PPC','wPLI'};
    band_index_map('Myogenic') = {'bilateral_autoregulation_efficiency', 'COx_CI_lower', 'COx_CI_upper', 'COx_stability', 'cardiac_cerebral_coupling', 'COx_single_hemisphere', 'COx_window_SE_rSO2', 'xcorr_max_corr', 'xcorr_optimal_lag', 'xcorr_asymmetry', 'xcorr_width','SE_power', 'rSO2_power', 'power_ratio', 'spectral_centroid_SE', 'spectral_centroid_rSO2', 'bandwidth_SE', 'bandwidth_rSO2', 'rolloff_SE', 'rolloff_rSO2','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality','Lx','myogenic_reactivity', 'smooth_muscle_tone_SE', 'smooth_muscle_tone_rSO2', 'myogenic_autoregulation_index', 'vascular_compliance', 'myogenic_frequency_stability','phase_diff','amplitude_corr','phase_amp_coupling', 'transfer_gain','transfer_phase','transfer_coherence','transfer_phase_lead','transfer_gain_variability','transfer_phase_variability','wavelet_coherence','wavelet_phase','PLV','PPC','wPLI'};
    band_index_map('Respiratory') = {'autonomic_modulation_index', 'cardiac_cerebral_coupling','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality','respiratory_coupling_strength', 'respiratory_modulation_SE', 'respiratory_modulation_rSO2', 'respiratory_phase_lag_sec', 'dominant_respiratory_frequency', 'respiratory_rate_variability'};
    band_index_map('Cardiac') = {'cardiac_coupling_strength', 'pulse_amplitude_SE', 'pulse_amplitude_rSO2', 'HRV_proxy', 'cardiac_phase_lag_sec', 'pulse_wave_velocity_proxy', 'dominant_cardiac_frequency', 'estimated_heart_rate','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality'};
    band_index_map('VLF') = {'autoregulation_effectiveness', 'autoregulation_time_constant', 'bilateral_autoregulation_efficiency', 'HAI_COx', 'COx_CI_lower', 'COx_CI_upper', 'COx_stability', 'cardiac_cerebral_coupling', 'COx_single_hemisphere', 'COx_window_SE_rSO2', 'xcorr_max_corr', 'xcorr_optimal_lag', 'xcorr_asymmetry', 'xcorr_width','SE_power', 'rSO2_power', 'power_ratio', 'spectral_centroid_SE', 'spectral_centroid_rSO2', 'bandwidth_SE', 'bandwidth_rSO2', 'rolloff_SE', 'rolloff_rSO2','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality','Lx','phase_diff','amplitude_corr','phase_amp_coupling', 'transfer_gain','transfer_phase','transfer_coherence','transfer_phase_lead','transfer_gain_variability','transfer_phase_variability','autoregulation_status','phase_margin','gain_margin_db','autoregulation_efficiency','wavelet_coherence','wavelet_phase','PLV','PPC','wPLI'};
    band_index_map('LF') = {'COx_single_hemisphere', 'COx_window_SE_rSO2', 'xcorr_max_corr', 'xcorr_optimal_lag', 'xcorr_asymmetry', 'xcorr_width','SE_power', 'rSO2_power', 'power_ratio', 'spectral_centroid_SE', 'spectral_centroid_rSO2', 'bandwidth_SE', 'bandwidth_rSO2', 'rolloff_SE', 'rolloff_rSO2','SE_to_rSO2_causality', 'rSO2_to_SE_causality', 'bidirectional_causality','Lx','phase_diff','amplitude_corr','phase_amp_coupling', 'transfer_gain','transfer_phase','transfer_coherence','transfer_phase_lead','transfer_gain_variability','transfer_phase_variability','autoregulation_status','phase_margin','gain_margin_db','autoregulation_efficiency','wavelet_coherence','wavelet_phase','PLV','PPC','wPLI'};
    
    % CSV header létrehozása
    csv_data = [];
    headers = {'PatientID'};
    
    % Minden sáv és index kombinációhoz
    applicable_bands = keys(band_index_map);
    for i = 1:length(applicable_bands)
        band_name = applicable_bands{i};
        indices_for_band = band_index_map(band_name);
        
        fprintf('  Checking band: %s\n', band_name);
        if isfield(specific_indices, band_name) && ~isempty(specific_indices.(band_name))
            fprintf('    Found data for %s: %d records\n', band_name, length(specific_indices.(band_name)));
            
            % Minden index ehhez a sávhoz
            for j = 1:length(indices_for_band)
                index_name = indices_for_band{j};
                if strcmp(index_name, 'autonomic_modulation_index')
                    headers{end+1} = sprintf('AutonomicModulation_%s', band_name);
                elseif strcmp(index_name, 'autoregulation_effectiveness')
                    headers{end+1} = sprintf('AutoregulationEffectiveness_%s', band_name);
                elseif strcmp(index_name, 'autoregulation_time_constant')
                    headers{end+1} = sprintf('AutoregulationTimeConstant_%s', band_name);
                elseif strcmp(index_name, 'bilateral_autoregulation_efficiency')
                    headers{end+1} = sprintf('BilateralEfficiency_%s', band_name);
                elseif strcmp(index_name, 'HAI_COx')
                    headers{end+1} = sprintf('HAI_COx_%s', band_name);
                elseif strcmp(index_name, 'COx_CI_lower')
                    headers{end+1} = sprintf('COx_CI_Lower_%s', band_name);
                elseif strcmp(index_name, 'COx_CI_upper')
                    headers{end+1} = sprintf('COx_CI_Upper_%s', band_name);
                elseif strcmp(index_name, 'COx_stability')
                    headers{end+1} = sprintf('COx_Stability_%s', band_name);
                elseif strcmp(index_name, 'cardiac_cerebral_coupling')
                    headers{end+1} = sprintf('CardiacCerebralCoupling_%s', band_name);
                elseif strcmp(index_name, 'cardiac_coupling_strength')
                    headers{end+1} = sprintf('CardiacCouplingStrength_%s', band_name);
                elseif strcmp(index_name, 'pulse_amplitude_SE')
                    headers{end+1} = sprintf('PulseAmplitudeSE_%s', band_name);
                elseif strcmp(index_name, 'pulse_amplitude_rSO2')
                    headers{end+1} = sprintf('PulseAmplitudeRSO2_%s', band_name);
                elseif strcmp(index_name, 'HRV_proxy')
                    headers{end+1} = sprintf('HRV_Proxy_%s', band_name);
                elseif strcmp(index_name, 'cardiac_phase_lag_sec')
                    headers{end+1} = sprintf('CardiacPhaseLag_%s', band_name);
                elseif strcmp(index_name, 'pulse_wave_velocity_proxy')
                    headers{end+1} = sprintf('PulseWaveVelocity_%s', band_name);
                elseif strcmp(index_name, 'dominant_cardiac_frequency')
                    headers{end+1} = sprintf('DominantCardiacFreq_%s', band_name);
                elseif strcmp(index_name, 'estimated_heart_rate')
                    headers{end+1} = sprintf('EstimatedHeartRate_%s', band_name);
                elseif strcmp(index_name, 'COx_single_hemisphere')
                    headers{end+1} = sprintf('COx_SingleHemisphere_%s', band_name);
                elseif strcmp(index_name, 'COx_window_SE_rSO2')
                    headers{end+1} = sprintf('COx_Window_SE_rSO2_%s', band_name);
                elseif strcmp(index_name, 'xcorr_max_corr')
                    headers{end+1} = sprintf('XCorr_MaxCorr_%s', band_name);
                elseif strcmp(index_name, 'xcorr_optimal_lag')
                    headers{end+1} = sprintf('XCorr_OptimalLag_%s', band_name);
                elseif strcmp(index_name, 'xcorr_asymmetry')
                    headers{end+1} = sprintf('XCorr_Asymmetry_%s', band_name);
                elseif strcmp(index_name, 'xcorr_width')
                    headers{end+1} = sprintf('XCorr_Width_%s', band_name);
               
                elseif strcmp(index_name, 'CVRI')
                    headers{end+1} = sprintf('CVRI');
                elseif strcmp(index_name, 'SE_alpha')
                    headers{end+1} = sprintf('SE_alpha');
                elseif strcmp(index_name, 'cross_DFA')
                    headers{end+1} = sprintf('cross_DFA');
                elseif strcmp(index_name, 'rSO2_alpha')
                    headers{end+1} = sprintf('rSO2_alpha');
                elseif strcmp(index_name, 'estimated_CCP')
                    headers{end+1} = sprintf('EstimatedCCP');
                elseif strcmp(index_name, 'estimated_RAP')
                    headers{end+1} = sprintf('EstimatedRAP');
                elseif strcmp(index_name, 'fractal_dim_SE')
                    headers{end+1} = sprintf('FractalDim_SE');
                elseif strcmp(index_name, 'fractal_dim_rSO2')
                    headers{end+1} = sprintf('FractalDim_rSO2');


                elseif strcmp(index_name, 'vasomotion_strength_SE')
                    headers{end+1} = sprintf('VasomotionStrength_SE');
                elseif strcmp(index_name, 'vasomotion_strength_rSO2')
                    headers{end+1} = sprintf('VasomotionStrength_rSO2');
                elseif strcmp(index_name, 'vasomotion_coupling')
                    headers{end+1} = sprintf('VasomotionCoupling');
                elseif strcmp(index_name, 'endothelial_dysfunction_index')
                    headers{end+1} = sprintf('EndothelialDysfunctionIndex');
                elseif strcmp(index_name, 'vasomotion_regularity')
                    headers{end+1} = sprintf('VasomotionRegularity');
                elseif strcmp(index_name, 'SE_power')
                    headers{end+1} = sprintf('SE_Power_%s', band_name);
                elseif strcmp(index_name, 'rSO2_power')
                    headers{end+1} = sprintf('rSO2_Power_%s', band_name);
                elseif strcmp(index_name, 'power_ratio')
                    headers{end+1} = sprintf('PowerRatio_%s', band_name);
                elseif strcmp(index_name, 'spectral_centroid_SE')
                    headers{end+1} = sprintf('SpectralCentroid_SE_%s', band_name);
                elseif strcmp(index_name, 'spectral_centroid_rSO2')
                    headers{end+1} = sprintf('SpectralCentroid_rSO2_%s', band_name);
                elseif strcmp(index_name, 'bandwidth_SE')
                    headers{end+1} = sprintf('Bandwidth_SE_%s', band_name);
                elseif strcmp(index_name, 'bandwidth_rSO2')
                    headers{end+1} = sprintf('Bandwidth_rSO2_%s', band_name);
                elseif strcmp(index_name, 'rolloff_SE')
                    headers{end+1} = sprintf('Rolloff_SE_%s', band_name);
                elseif strcmp(index_name, 'rolloff_rSO2')
                    headers{end+1} = sprintf('Rolloff_rSO2_%s', band_name);

                elseif strcmp(index_name, 'SE_to_rSO2_causality')
                    headers{end+1} = sprintf('SE_to_rSO2_Causality_%s', band_name);
                elseif strcmp(index_name, 'rSO2_to_SE_causality')
                    headers{end+1} = sprintf('rSO2_to_SE_Causality_%s', band_name);
                elseif strcmp(index_name, 'bidirectional_causality')
                    headers{end+1} = sprintf('Bidirectional_Causality_%s', band_name);

                elseif strcmp(index_name, 'HRx')
                    headers{end+1} = sprintf('HRx');
                elseif strcmp(index_name, 'TOHRx')
                    headers{end+1} = sprintf('TOHRx');

                elseif strcmp(index_name, 'hurst_SE')
                    headers{end+1} = sprintf('Hurst_SE');
                elseif strcmp(index_name, 'hurst_rSO2')
                    headers{end+1} = sprintf('Hurst_rSO2');

               elseif strcmp(index_name, 'Lx')
                    headers{end+1} = sprintf('Lx_%s', band_name);

               elseif strcmp(index_name, 'SE_MSE')
                   headers{end+1} = sprintf('SE_MSE');
               elseif strcmp(index_name, 'rSO2_MSE')
                   headers{end+1} = sprintf('rSO2_MSE');
               elseif strcmp(index_name, 'complexity_index')
                   headers{end+1} = sprintf('ComplexityIndex');
               elseif strcmp(index_name, 'mutual_information')
                    headers{end+1} = sprintf('MutualInformation');
               elseif strcmp(index_name, 'Mx')
                    headers{end+1} = sprintf('Mx');

                elseif strcmp(index_name, 'myogenic_reactivity')
                    headers{end+1} = sprintf('MyogenicReactivity');
                elseif strcmp(index_name, 'smooth_muscle_tone_SE')
                    headers{end+1} = sprintf('SmoothMuscleTone_SE');
                elseif strcmp(index_name, 'smooth_muscle_tone_rSO2')
                    headers{end+1} = sprintf('SmoothMuscleTone_rSO2');
                elseif strcmp(index_name, 'myogenic_autoregulation_index')
                    headers{end+1} = sprintf('MyogenicAutoregulationIndex');
                
                elseif strcmp(index_name, 'myogenic_frequency_stability')
                    headers{end+1} = sprintf('MyogenicFrequencyStability');

                elseif strcmp(index_name, 'sympathetic_tone_SE')
                    headers{end+1} = sprintf('SympatheticTone_SE');
                elseif strcmp(index_name, 'sympathetic_tone_rSO2')
                    headers{end+1} = sprintf('SympatheticTone_rSO2');
                elseif strcmp(index_name, 'autonomic_balance')
                    headers{end+1} = sprintf('AutonomicBalance');
                elseif strcmp(index_name, 'neurogenic_coupling_efficiency')
                    headers{end+1} = sprintf('NeurogenicCouplingEfficiency');
                elseif strcmp(index_name, 'neurogenic_burst_frequency')
                    headers{end+1} = sprintf('NeurogenicBurstFrequency');
                elseif strcmp(index_name, 'neurogenic_burst_duration')
                    headers{end+1} = sprintf('NeurogenicBurstDuration');
                elseif strcmp(index_name, 'normalized_mutual_information')
                    headers{end+1} = sprintf('NormalizedMutualInformation');

                elseif strcmp(index_name, 'COx_temporal_variability')
                    headers{end+1} = sprintf('COx_TemporalVariability');
                elseif strcmp(index_name, 'PRx')
                    headers{end+1} = sprintf('PRx');
                elseif strcmp(index_name, 'autoregulation_range')
                    headers{end+1} = sprintf('AutoregulationRange');
                elseif strcmp(index_name, 'FRx')
                    headers{end+1} = sprintf('FRx');
                elseif strcmp(index_name, 'vascular_compliance')
                    if strcmp(band_name, 'Endothelial')
                        headers{end+1} = sprintf('VascularCompliance_Global');
                    elseif strcmp(band_name, 'Myogenic') 
                        headers{end+1} = sprintf('VascularCompliance_Myogenic');
                    else
                        headers{end+1} = sprintf('VascularCompliance_%s', band_name);
                    end
                
                elseif strcmp(index_name, 'respiratory_coupling_strength')
                    headers{end+1} = sprintf('RespiratoryCouplingStrength');
                elseif strcmp(index_name, 'respiratory_modulation_SE')
                    headers{end+1} = sprintf('RespiratoryModulation_SE');
                elseif strcmp(index_name, 'respiratory_modulation_rSO2')
                    headers{end+1} = sprintf('RespiratoryModulation_rSO2');
                elseif strcmp(index_name, 'respiratory_phase_lag_sec')
                    headers{end+1} = sprintf('RespiratoryPhaseLag');
                elseif strcmp(index_name, 'dominant_respiratory_frequency')
                    headers{end+1} = sprintf('DominantRespiratoryFreq');
                elseif strcmp(index_name, 'respiratory_rate_variability')
                    headers{end+1} = sprintf('RespiratoryRateVariability');

                elseif strcmp(index_name, 'surrogate_p_value')
                    headers{end+1} = sprintf('SurrogatePValue');
                elseif strcmp(index_name, 'surrogate_significance')
                    headers{end+1} = sprintf('SurrogateSignificance');

                elseif strcmp(index_name, 'Sx')
                    headers{end+1} = sprintf('Sx');

                elseif strcmp(index_name, 'phase_diff')
                    headers{end+1} = sprintf('PhaseDiff_%s', band_name);
                elseif strcmp(index_name, 'amplitude_corr')
                    headers{end+1} = sprintf('AmplitudeCorr_%s', band_name);
                elseif strcmp(index_name, 'phase_amp_coupling')
                    headers{end+1} = sprintf('PhaseAmpCoupling_%s', band_name);
                 elseif strcmp(index_name, 'transfer_gain')
                    headers{end+1} = sprintf('TransferGain_%s', band_name);
                elseif strcmp(index_name, 'transfer_phase')
                    headers{end+1} = sprintf('TransferPhase_%s', band_name);
                elseif strcmp(index_name, 'transfer_coherence')
                    headers{end+1} = sprintf('TransferCoherence_%s', band_name);
                elseif strcmp(index_name, 'transfer_phase_lead')
                    headers{end+1} = sprintf('TransferPhaseLead_%s', band_name);
                elseif strcmp(index_name, 'transfer_gain_variability')
                    headers{end+1} = sprintf('TransferGainVariability_%s', band_name);
                elseif strcmp(index_name, 'transfer_phase_variability')
                    headers{end+1} = sprintf('TransferPhaseVariability_%s', band_name);
                elseif strcmp(index_name, 'autoregulation_status')
                    headers{end+1} = sprintf('AutoregulationStatus_%s', band_name);
                elseif strcmp(index_name, 'phase_margin')
                    headers{end+1} = sprintf('PhaseMargin_%s', band_name);
                elseif strcmp(index_name, 'gain_margin_db')
                    headers{end+1} = sprintf('GainMargin_%s', band_name);
                elseif strcmp(index_name, 'autoregulation_efficiency')
                    headers{end+1} = sprintf('AutoregulationEfficiency_%s', band_name);
                elseif strcmp(index_name, 'wavelet_coherence')
                    headers{end+1} = sprintf('WaveletCoherence_%s', band_name);
                elseif strcmp(index_name, 'wavelet_phase')
                    headers{end+1} = sprintf('WaveletPhase_%s', band_name);
                elseif strcmp(index_name, 'PLV')
                    headers{end+1} = sprintf('PLV_%s', band_name);
                elseif strcmp(index_name, 'PPC')
                    headers{end+1} = sprintf('PPC_%s', band_name);
                elseif strcmp(index_name, 'wPLI')
                    headers{end+1} = sprintf('wPLI_%s', band_name);
                end
            end
        else
            fprintf('    No data found for %s\n', band_name);
        end
    end
    
    % Ha nincsenek alkalmazható indexek
    if length(headers) == 1
        fprintf('HIBA: Nem találhatók indexek!\n');
        fprintf('Elérhető sávok: %s\n', strjoin(fieldnames(specific_indices), ', '));
        return;
    end
    
    % Minden páciens adatainak összeállítása
    for p = 1:length(all_patients)
        patient_id = all_patients(p);
        patient_row = {patient_id};
        
        % Minden sáv és index kombinációhoz
        for i = 1:length(applicable_bands)
            band_name = applicable_bands{i};
            indices_for_band = band_index_map(band_name);
            
            if isfield(specific_indices, band_name) && ~isempty(specific_indices.(band_name))
                % Keressük meg ezt a pácienst ebben a sávban
                band_indices = specific_indices.(band_name);
                patient_indices = band_indices([band_indices.PatientID] == patient_id);
                
                if ~isempty(patient_indices)
                    % Minden index értékét hozzáadjuk
                    for j = 1:length(indices_for_band)
                        index_name = indices_for_band{j};
                        if isfield(patient_indices, index_name)
                            patient_row{end+1} = patient_indices.(index_name);
                        else
                            patient_row{end+1} = NaN;
                        end
                    end
                else
                    % Páciens nem volt ebben a sávban - minden indexhez NaN
                    for j = 1:length(indices_for_band)
                        patient_row{end+1} = NaN;
                    end
                end
            end
        end
        
        csv_data = [csv_data; patient_row];
    end
    
    % CSV létrehozása és mentése
    if ~isempty(csv_data) && size(csv_data, 2) > 1
        csv_table = cell2table(csv_data, 'VariableNames', headers);
        writetable(csv_table, 'specific_indices_matrix.csv');
        fprintf('Specifikus indexek exportálva: specific_indices_matrix.csv\n');
        fprintf('Mátrix méret: %d páciensek × %d indexek\n', height(csv_table), width(csv_table)-1);
        
        % Összefoglaló statisztikák
        fprintf('\n=== INDEXEK ÖSSZEFOGLALÓJA ===\n');
        for i = 2:length(headers) % Skip PatientID
            index_name = headers{i};
            values = [csv_data{:,i}];
            valid_values = values(~isnan(values));
            
            fprintf('%s:\n', index_name);
            fprintf('  Valid páciensek: %d/%d (%.1f%%)\n', length(valid_values), length(values), length(valid_values)/length(values)*100);
            if ~isempty(valid_values)
                fprintf('  Átlag: %.4f ± %.4f\n', mean(valid_values), std(valid_values));
                fprintf('  Tartomány: %.4f - %.4f\n', min(valid_values), max(valid_values));
            end
            fprintf('\n');
        end
    else
        fprintf('Nincs exportálható index adat vagy csak PatientID oszlop van.\n');
    end
end


function export_specific_indices_to_csv2(specific_indices, all_patients)
    % CSV export: Csak a kiválasztott metrikák
    % Formátum: sorok = metrikák, oszlopok = páciensek
    
    % Kiválasztott metrikák definíciója
    selected_metrics = {
        {'COx_single_hemisphere', 'Endothelial', 'COx_SingleHemisphere_Endothelial'};
        {'vasomotion_strength_SE', 'Endothelial', 'VasomotionStrength_SE'};
        {'rSO2_power', 'Endothelial', 'rSO2_Power_Endothelial'};
        {'power_ratio', 'Endothelial', 'PowerRatio_Endothelial'};
        {'Lx', 'Endothelial', 'Lx_Endothelial'};
        {'transfer_gain_variability', 'LF', 'TransferGainVariability_LF'};
        {'transfer_phase_variability', 'Myogenic', 'TransferPhaseVariability_Myogenic'};
        {'autonomic_modulation_index', 'Neurogenic', 'AutonomicModulation_Neurogenic'};
        {'power_ratio', 'Neurogenic', 'PowerRatio_Neurogenic'};
        {'Lx', 'Neurogenic', 'Lx_Neurogenic'};
        {'amplitude_corr', 'Neurogenic', 'AmplitudeCorr_Neurogenic'};
        {'wavelet_coherence', 'Neurogenic', 'WaveletCoherence_Neurogenic'};
        {'rSO2_to_SE_causality', 'Respiratory', 'rSO2_to_SE_Causality_Respiratory'};
        {'bidirectional_causality', 'Respiratory', 'Bidirectional_Causality_Respiratory'};
        {'autoregulation_time_constant', 'VLF', 'AutoregulationTimeConstant_VLF'};
        {'COx_single_hemisphere', 'VLF', 'COx_SingleHemisphere_VLF'};
        {'power_ratio', 'VLF', 'PowerRatio_VLF'};
        {'transfer_phase_variability', 'VLF', 'TransferPhaseVariability_VLF'};
    };
    
    % Inicializálás
    n_metrics = length(selected_metrics);
    n_patients = length(all_patients);
    
    % Adatmátrix létrehozása (sorok = metrikák, oszlopok = páciensek)
    data_matrix = NaN(n_metrics, n_patients);
    metric_names = cell(n_metrics, 1);
    
    % Metrikák kinyerése
    for m = 1:n_metrics
        metric_field = selected_metrics{m}{1};
        band_name = selected_metrics{m}{2};
        metric_label = selected_metrics{m}{3};
        metric_names{m} = metric_label;
        
        % Ellenőrizzük, hogy létezik-e ez a sáv az adatokban
        if isfield(specific_indices, band_name) && ~isempty(specific_indices.(band_name))
            band_data = specific_indices.(band_name);
            
            % Minden pácienshez keressük meg az értéket
            for p = 1:n_patients
                patient_id = all_patients(p);
                
                % Keressük meg a pácienst a sáv adataiban
                patient_idx = find([band_data.PatientID] == patient_id);
                
                if ~isempty(patient_idx) && isfield(band_data(patient_idx), metric_field)
                    data_matrix(m, p) = band_data(patient_idx).(metric_field);
                end
            end
        end
    end
    
    % CSV tábla létrehozása
    % Első oszlop: metrika nevek
    csv_table = table(metric_names, 'VariableNames', {'Metric'});
    
    % Páciens oszlopok hozzáadása
    for p = 1:n_patients
        patient_col_name = sprintf('Patient_%d', all_patients(p));
        csv_table.(patient_col_name) = data_matrix(:, p);
    end
    
    % CSV mentése
    writetable(csv_table, 'selected_metrics_export.csv');
    
    fprintf('\n=== KIVÁLASZTOTT METRIKÁK EXPORTÁLVA ===\n');
    fprintf('Fájl: selected_metrics_export.csv\n');
    fprintf('Mátrix méret: %d metrika × %d páciens\n', n_metrics, n_patients);
    
    % Összefoglaló statisztikák
    fprintf('\n=== METRIKA STATISZTIKÁK ===\n');
    for m = 1:n_metrics
        metric_values = data_matrix(m, :);
        valid_values = metric_values(~isnan(metric_values));
        
        fprintf('\n%s:\n', metric_names{m});
        fprintf('  Valid páciensek: %d/%d (%.1f%%)\n', ...
            length(valid_values), n_patients, length(valid_values)/n_patients*100);
        
        if ~isempty(valid_values)
            fprintf('  Átlag: %.4f ± %.4f\n', mean(valid_values), std(valid_values));
            fprintf('  Medián: %.4f\n', median(valid_values));
            fprintf('  Tartomány: [%.4f, %.4f]\n', min(valid_values), max(valid_values));
        end
    end
    
    % Transzponált verzió mentése (ha szükséges)
    transposed_data = array2table(data_matrix', 'VariableNames', metric_names);
    transposed_data.PatientID = all_patients;
    % PatientID-t az első oszlopba tesszük
    transposed_data = transposed_data(:, [end, 1:end-1]);
    
    writetable(transposed_data, 'selected_metrics_transposed.csv');
    fprintf('\n=== TRANSZPONÁLT VERZIÓ IS MENTVE ===\n');
    fprintf('Fájl: selected_metrics_transposed.csv\n');
    fprintf('Mátrix méret: %d páciens × %d metrika\n', n_patients, n_metrics);
end
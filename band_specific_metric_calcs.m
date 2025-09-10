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

% Debug információ - MOST JÓ HELYEN!
fprintf('\n=== ADAT DEBUG ===\n');
fprintf('SE oszlop típus: %s\n', class(data.SE));
fprintf('SE tartomány: %.1f - %.1f\n', nanmin(data.SE), nanmax(data.SE));
fprintf('SE NaN arány: %.1f%%\n', sum(isnan(data.SE))/length(data.SE)*100);
fprintf('MAP tartomány: %.1f - %.1f\n', nanmin(data.MAP), nanmax(data.MAP));
fprintf('rSO2 tartomány: %.1f - %.1f\n', nanmin(data.oper_side_oxig), nanmax(data.oper_side_oxig));

%% 2. SÁVONKÉNTI ADATFELKÉSZÍTÉS
processed_data = struct();

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
        
        % Valid arány ellenőrzés
        se_valid_ratio = sum(~isnan(SE_signal)) / length(SE_signal);
        rso2_valid_ratio = sum(~isnan(rSO2_signal)) / length(rSO2_signal);
        map_valid_ratio = sum(~isnan(MAP_signal)) / length(MAP_signal);
        
        if se_valid_ratio < min_valid_ratio || rso2_valid_ratio < min_valid_ratio || map_valid_ratio < min_valid_ratio
            fprintf('  Páciens %d: elégtelen adatminőség (SE: %.1f%%, rSO2: %.1f%%, MAP: %.1f%%)\n', ...
                    patient_id, se_valid_ratio*100, rso2_valid_ratio*100, map_valid_ratio*100);
            continue;
        end
        
        % 2. Gap filling
        max_gap_points = round(band_params.max_gap_seconds * band_params.target_fs);
        SE_cleaned = smart_gap_filling(SE_signal, max_gap_points);
        rSO2_cleaned = smart_gap_filling(rSO2_signal, max_gap_points);
        MAP_cleaned = smart_gap_filling(MAP_signal, max_gap_points);
        
        try
            % 3. Sáv-specifikus interpoláció
            patient_data_for_interp = patient_data;
            patient_data_for_interp.SE = SE_cleaned;
            patient_data_for_interp.oper_side_oxig = rSO2_cleaned;
            patient_data_for_interp.MAP = MAP_cleaned;
            
            interpolated_data = interpolate_data(patient_data_for_interp, ...
                                                band_params.n_points, ...
                                                band_params.optimal_duration);
            
            SE_interp = interpolated_data.SE;
            rSO2_interp = interpolated_data.rSO2;
            MAP_interp = interpolated_data.MAP;
            
            % 4. Missing values kezelése
            SE_clean = fillmissing(SE_interp, 'nearest');
            rSO2_clean = fillmissing(rSO2_interp, 'nearest');
            MAP_clean = fillmissing(MAP_interp, 'nearest');
            
            % 4.5. FIZIOLÓGIAI HATÁROK ELLENŐRZÉSE (RAW ADATOKON!)
            physio_bounds = band_params.physio_bounds;
            SE_physio_valid = SE_clean >= physio_bounds(1,1) & SE_clean <= physio_bounds(1,2);
            rSO2_physio_valid = rSO2_clean >= physio_bounds(2,1) & rSO2_clean <= physio_bounds(2,2);
            MAP_physio_valid = MAP_clean >= physio_bounds(3,1) & MAP_clean <= physio_bounds(3,2);
            
            % Valid arány RAW adatokon
            physio_valid_ratio = sum(SE_physio_valid & rSO2_physio_valid & MAP_physio_valid) / length(SE_clean);
            
            if physio_valid_ratio < min_valid_ratio
                fprintf('  Páciens %d: fiziológiai határok sértése (%.1f%% valid)\n', ...
                        patient_id, physio_valid_ratio*100);
                continue;
            end
            
            % 5. Detrending (sáv-specifikus) - CSAK ezután!
            SE_clean = advanced_detrend(SE_clean, band_params.detrend_method, 3);
            rSO2_clean = advanced_detrend(rSO2_clean, band_params.detrend_method, 3);
            MAP_clean = advanced_detrend(MAP_clean, band_params.detrend_method, 3);
            
            % 6. Sáv-specifikus high-pass szűrés (egyszerűsített)
            cutoff_freq = max(0.001, band_params.range(1) * 0.5);
            if length(SE_clean) > 10 && cutoff_freq < band_params.target_fs/2
                % Egyszerű Butterworth high-pass szűrő (gyorsabb, stabilabb)
                try
                    [b, a] = butter(1, cutoff_freq/(band_params.target_fs/2), 'high');
                    SE_clean = filtfilt(b, a, SE_clean);
                    rSO2_clean = filtfilt(b, a, rSO2_clean);
                    MAP_clean = filtfilt(b, a, MAP_clean);
                catch
                    fprintf('    Szűrés kihagyva (nem kritikus)\n');
                end
            end
            
            % 7. Final validity check (feldolgozott adatok)
            % Most már csak finite értékeket ellenőrzünk
            SE_valid = isfinite(SE_clean);
            rSO2_valid = isfinite(rSO2_clean);
            MAP_valid = isfinite(MAP_clean);
            
            % Valid pontok aránya a tisztított adatban
            valid_ratio = sum(SE_valid & rSO2_valid & MAP_valid) / length(SE_clean);
            
            % Debug: értékek ellenőrzése
            if p <= 3  % Első 3 páciensre debug info
                fprintf('    DEBUG - SE: %.1f-%.1f, rSO2: %.1f-%.1f, MAP: %.1f-%.1f\n', ...
                        nanmin(SE_clean), nanmax(SE_clean), nanmin(rSO2_clean), nanmax(rSO2_clean), ...
                        nanmin(MAP_clean), nanmax(MAP_clean));
            end
            
            if valid_ratio >= min_valid_ratio
                % Készítjük fel az adatokat a metrika számításokhoz
                processed_patient = struct();
                processed_patient.PatientID = patient_id;
                processed_patient.SE_clean = SE_clean;
                processed_patient.rSO2_clean = rSO2_clean;
                processed_patient.MAP_clean = MAP_clean;
                processed_patient.valid_ratio = valid_ratio;
                processed_patient.sampling_rate = band_params.target_fs;
                processed_patient.freq_range = band_params.range;
                processed_patient.data_length_sec = length(SE_clean) / band_params.target_fs;
                
                band_data.patient_data = [band_data.patient_data; processed_patient];
                band_data.valid_patients = [band_data.valid_patients; patient_id];
                
                fprintf('  Páciens %d: siker (%.1f%% valid, %.1f sec)\n', ...
                        patient_id, valid_ratio*100, processed_patient.data_length_sec);
            else
                fprintf('  Páciens %d: fiziológiai határok sértése (%.1f%% valid)\n', ...
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
        for i = 1:length(band_data.patient_data)
            all_SE = [all_SE; band_data.patient_data(i).SE_clean];
            all_rSO2 = [all_rSO2; band_data.patient_data(i).rSO2_clean];
            all_MAP = [all_MAP; band_data.patient_data(i).MAP_clean];
        end
        
        fprintf('  SE átlag: %.1f ± %.1f mmHg\n', mean(all_SE), std(all_SE));
        fprintf('  rSO2 átlag: %.1f ± %.1f %%\n', mean(all_rSO2), std(all_rSO2));
        fprintf('  MAP átlag: %.1f ± %.1f mmHg\n', mean(all_MAP), std(all_MAP));
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
fprintf('Következő lépés: metrika számítások hozzáadása.\n');

%% ========================================================================
%% SEGÉDFÜGGVÉNYEK HELYE
%% ========================================================================
% A következő függvényeknek a functions/ vagy functions/helper/ mappában kell lenniük:
% - smart_gap_filling.m (functions/helper/ mappa)
% - interpolate_data.m (functions/ mappa - MAP támogatással frissítve)
% - advanced_detrend.m (functions/ mappa)
%
% FIXEK ALKALMAZVA:
% - nanmin/nanmax használata min(...,'omitnan') helyett (régi MATLAB kompatibilitás)
% - Fiziológiai határok lazítva: SE: 0-200, rSO2: 0-100, MAP: 30-200
% - High-pass szűrő lecserélve egyszerű Butterworth-ra (gyorsabb)
% - Debug információ hozzáadva az első páciensekhez
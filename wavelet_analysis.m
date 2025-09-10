% BAND_SPECIFIC_WAVELET_ANALYSIS - Sávonként optimalizált analízis
% SE és rSO2 közötti wavelet coupling analízis frekvencia sávonként

clear; clc;

%% Paraméterek
data_file = 'df_unnormalized.csv';
target_oppart = 4;  % Clamp fázis
min_valid_ratio = 0.70;  % Min 70% valid adat kell

fprintf('=== SÁV-SPECIFIKUS WAVELET ANALÍZIS ===\n');

%% Frekvencia sávok definíciója optimális paraméterekkel

freq_bands = struct();

% Endothelial - lassú, hosszú mérés
freq_bands.Endothelial = struct( ...
    'range', [0.003, 0.02], ...
    'optimal_duration', 1000, ...       % 1000 sec (16 perc) minimum
    'target_fs', 0.1, ...               % 0.1 Hz sampling (10 sec/pont)
    'n_points', 100, ...                % 100 pont = 1000 sec
    'min_cycles', 3, ...                % Min 3 teljes ciklus
    'max_gap_seconds', 50, ...          % 50 sec gap tolerancia (lassú folyamat)
    'detrend_method', 'emd', ...        % EMD detrending (non-linear trends)
    'detrend_order', [], ...            % EMD-hez nincs order
    'artifact_z_threshold', 3.0, ...    % Enyhébb artifact detection
    'artifact_gradient_threshold', 0.5, ... % Lassú változások OK
    'physio_bounds', [30, 85; 40, 95]); % [SE_min, SE_max; rSO2_min, rSO2_max]

% Neurogenic - közepes
freq_bands.Neurogenic = struct( ...
    'range', [0.02, 0.06], ...
    'optimal_duration', 500, ...        % 500 sec (8 perc)
    'target_fs', 0.2, ...               % 0.2 Hz sampling (5 sec/pont)
    'n_points', 100, ...
    'min_cycles', 4, ...
    'max_gap_seconds', 25, ...          % 25 sec gap tolerancia
    'detrend_method', 'emd', ... % Polynomial detrending
    'detrend_order', 3, ...             % 3rd order polynomial
    'artifact_z_threshold', 3.2, ...    
    'artifact_gradient_threshold', 0.8, ...
    'physio_bounds', [30, 85; 40, 95]);

% Myogenic - gyors, rövid mérés is elég
freq_bands.Myogenic = struct( ...
    'range', [0.06, 0.15], ...
    'optimal_duration', 300, ...        % 300 sec (5 perc)
    'target_fs', 0.5, ...               % 0.5 Hz sampling (2 sec/pont)
    'n_points', 150, ...
    'min_cycles', 5, ...
    'max_gap_seconds', 12, ...          % 12 sec gap tolerancia
    'detrend_method', 'emd', ... % Quadratic detrending
    'detrend_order', 2, ...             % 2nd order polynomial
    'artifact_z_threshold', 3.5, ...    
    'artifact_gradient_threshold', 1.2, ...
    'physio_bounds', [30, 85; 40, 95]);

% Respiratory - nagyon gyors
freq_bands.Respiratory = struct( ...
    'range', [0.15, 0.4], ...
    'optimal_duration', 200, ...        % 200 sec (3 perc)
    'target_fs', 1.0, ...               % 1 Hz sampling (1 sec/pont)
    'n_points', 200, ...
    'min_cycles', 6, ...
    'max_gap_seconds', 6, ...           % 6 sec gap tolerancia
    'detrend_method', 'linear', ...     % Linear detrending (gyors változások)
    'detrend_order', 1, ...             % Linear
    'artifact_z_threshold', 4.0, ...    % Szigorúbb artifact detection
    'artifact_gradient_threshold', 2.0, ...
    'physio_bounds', [30, 85; 40, 95]);

% Cardiac - extrém gyors (ha elérhető)
freq_bands.Cardiac = struct( ...
    'range', [0.4, 2.0], ...
    'optimal_duration', 120, ...        % 120 sec (2 perc)
    'target_fs', 5.0, ...               % 5 Hz sampling
    'n_points', 600, ...
    'min_cycles', 10, ...
    'max_gap_seconds', 3, ...           % 3 sec gap tolerancia (nagyon gyors)
    'detrend_method', 'linear', ...     % Linear detrending
    'detrend_order', 1, ...             
    'artifact_z_threshold', 4.5, ...    % Legszűrőbb artifact detection
    'artifact_gradient_threshold', 3.0, ...
    'physio_bounds', [30, 85; 40, 95]);

band_names = fieldnames(freq_bands);
%% 1. Adatok betöltése (közös rész)
data = readtable(data_file);
% Próbáld mindkét oszlopnevet:
if ismember('oppart', data.Properties.VariableNames)
    phase_data = data(data.oppart == target_oppart, :);
elseif ismember('oper_phase', data.Properties.VariableNames)
    phase_data = data(data.oper_phase == target_oppart, :);
else
    % Debug: oszlopnevek kiírása
    fprintf('Elérhető oszlopok:\n');
    disp(data.Properties.VariableNames);
    error('Nem található oppart vagy oper_phase oszlop!');
end
patients = unique(phase_data.Identifier);
n_patients = length(patients);

fprintf('Eredeti páciensek száma: %d\n', n_patients);

%% 2. Sávonkénti analízis
results_all_bands = struct();
patient_results = table();

for band_idx = 1:length(band_names)
    band_name = band_names{band_idx};
    band_params = freq_bands.(band_name);
    
    fprintf('\n=== %s SÁV ANALÍZIS ===\n', band_name);
    fprintf('Frekvencia: %.3f-%.3f Hz\n', band_params.range(1), band_params.range(2));
    fprintf('Target fs: %.3f Hz, Pontok: %d\n', band_params.target_fs, band_params.n_points);
    
    % Sáv-specifikus eredmények tárolása
    band_results = struct();
    band_results.coherence_values = [];
    band_results.PLV_values = [];
    band_results.CMP_values = [];
    band_results.COx_values = [];
    band_results.valid_patients = [];
    
    %% Páciens loop ezen sávra
    for p = 1:n_patients
        patient_id = patients(p);
        patient_data = phase_data(phase_data.Identifier == patient_id, :);
        
        % 1. Adatminőség ellenőrzés
        SE_signal = patient_data.MAP;
        rSO2_signal = patient_data.oper_side_oxig;
        
        % Valid arány
        se_valid_ratio = sum(~isnan(SE_signal)) / length(SE_signal);
        rso2_valid_ratio = sum(~isnan(rSO2_signal)) / length(rSO2_signal);
        
        if se_valid_ratio < min_valid_ratio || rso2_valid_ratio < min_valid_ratio
            continue;
        end
        
        % 2. Gap filling
        max_gap_points = round(band_params.max_gap_seconds * band_params.target_fs);
        SE_cleaned = smart_gap_filling(SE_signal, max_gap_points);
        rSO2_cleaned = smart_gap_filling(rSO2_signal, max_gap_points);
        patient_data_cleaned = patient_data;
        patient_data_cleaned.SE = SE_cleaned;
        patient_data_cleaned.oper_side_oxig = rSO2_cleaned;
        % 3. Sáv-specifikus interpoláció
        try
            % Egyedi interpoláció erre a sávra
            interpolated_data = interpolate_data(patient_data_cleaned, ...
                                                   band_params.n_points, ...
                                                   band_params.optimal_duration);
            
            SE_interp = interpolated_data.SE;
            rSO2_interp = interpolated_data.rSO2;
            
            % 4. Preprocessing
            SE_clean = fillmissing(SE_interp, 'nearest');
            rSO2_clean = fillmissing(rSO2_interp, 'nearest');
            
            SE_clean = advanced_detrend(SE_clean, band_params.detrend_method, band_params.detrend_order);
            rSO2_clean = advanced_detrend(rSO2_clean, band_params.detrend_method, band_params.detrend_order);
            % Sáv-specifikus high-pass szűrés
            cutoff_freq = max(0.001, band_params.range(1) * 0.5);
            if length(SE_clean) > 10
                SE_clean = highpass(SE_clean, cutoff_freq, band_params.target_fs);
                rSO2_clean = highpass(rSO2_clean, cutoff_freq, band_params.target_fs);
            end
            
            % 5. Wavelet analízis (sáv-specifikus)
            [wcoh, wcs, f] = wcoherence(SE_clean, rSO2_clean, band_params.target_fs);
            
            % Frekvencia maszk erre a sávra
            freq_mask = f >= band_params.range(1) & f <= band_params.range(2);
            
            if sum(freq_mask) >= band_params.min_cycles
                % Sáv-specifikus metrikák számítása
                coh_in_band = wcoh(freq_mask);
                phase_in_band = angle(wcs(freq_mask));
                
                % Validálás és tisztítás
                coh_in_band = real(coh_in_band);
                coh_in_band = coh_in_band(isfinite(coh_in_band) & coh_in_band >= 0 & coh_in_band <= 1);
                phase_in_band = phase_in_band(isfinite(phase_in_band));
                
                if length(coh_in_band) >= 3
                    % Eredmények tárolása
                    median_coh = median(coh_in_band);
                    PLV_val = abs(mean(exp(1i * phase_in_band)));
                    CMP_val = angle(mean(exp(1i * phase_in_band)));
                    
                    band_results.coherence_values = [band_results.coherence_values; median_coh];
                    band_results.PLV_values = [band_results.PLV_values; PLV_val];
                    band_results.CMP_values = [band_results.CMP_values; CMP_val];
                    band_results.valid_patients = [band_results.valid_patients; patient_id];
                    

                    
                    % COx számítás
                    try
                        COx_vals = compute_COx_index(SE_clean, rSO2_clean);
                        if ~isempty(COx_vals)
                            band_results.COx_values = [band_results.COx_values; mean(COx_vals)];
                        end
                    catch
                        band_results.COx_values = [band_results.COx_values; NaN];
                    end

                    new_patient_row = table(patient_id, {band_name}, median_coh, PLV_val, ...
                                           rad2deg(CMP_val), mean(COx_vals), ...
                                           'VariableNames', {'PatientID', 'Band', 'Coherence', 'PLV', 'CMP_deg', 'COx'});
                    patient_results = [patient_results; new_patient_row];

                end
            end
            
        catch ME
            % Hibás páciens kihagyása
            continue;
        end
        
        % Progress
        if mod(length(band_results.valid_patients), 10) == 0
            fprintf('  %s: %d páciens elemezve...\n', band_name, length(band_results.valid_patients));
        end
    end
    writetable(patient_results, 'patient_wavelet_results.csv');
    fprintf('Betegenkénti eredmények: patient_wavelet_results.csv\n');
    
    %% Sáv összegzés
    n_valid = length(band_results.valid_patients);
    fprintf('%s sáv: %d/%d páciens sikerült (%.1f%%)\n', ...
            band_name, n_valid, n_patients, n_valid/n_patients*100);
    
    if n_valid >= 5  % Minimum 5 páciens kell
        % Statisztikák
        band_results.mean_coherence = mean(band_results.coherence_values);
        band_results.std_coherence = std(band_results.coherence_values);
        band_results.mean_PLV = mean(band_results.PLV_values);
        band_results.mean_CMP = angle(mean(exp(1i * band_results.CMP_values)));
        band_results.mean_COx = mean(band_results.COx_values(~isnan(band_results.COx_values)));
        
        fprintf('  Átlag Coherence: %.3f ± %.3f\n', band_results.mean_coherence, band_results.std_coherence);
        fprintf('  Átlag PLV: %.3f\n', band_results.mean_PLV);
        fprintf('  Átlag CMP: %.1f°\n', rad2deg(band_results.mean_CMP));
        fprintf('  Átlag COx: %.3f\n', band_results.mean_COx);
    else
        fprintf('  ELÉGTELENÜL ADAT - sáv kihagyva\n');
        band_results = [];
    end
    
    % Eredmények mentése
    results_all_bands.(band_name) = band_results;
    
    %% Sáv-specifikus plot
    if ~isempty(band_results) && n_valid >= 5
        create_band_specific_plot(band_results, band_name, band_params);
    end
end

%% 3. Összesített eredmények és összehasonlítás
fprintf('\n=== ÖSSZESÍTETT EREDMÉNYEK ===\n');

% Összehasonlító táblázat
summary_table = table();
for band_idx = 1:length(band_names)
    band_name = band_names{band_idx};
    if isfield(results_all_bands, band_name) && ~isempty(results_all_bands.(band_name))
        result = results_all_bands.(band_name);
        
        new_row = table({band_name}, result.mean_coherence, result.mean_PLV, ...
                       rad2deg(result.mean_CMP), result.mean_COx, ...
                       length(result.valid_patients), ...
                       'VariableNames', {'Band', 'Coherence', 'PLV', 'CMP_deg', 'COx', 'N_patients'});
        summary_table = [summary_table; new_row];
    end
end

disp(summary_table);

% Összehasonlító plot
create_summary_comparison_plot(results_all_bands, band_names);

% CSV export
writetable(summary_table, 'band_specific_results.csv');
plot_wavelet_results('band_specific_results.csv');
fprintf('\nSáv-specifikus analízis befejezve!\n');
fprintf('Eredmények: band_specific_results.csv\n');



function signal_detrended = advanced_detrend(signal, method, order)
% ADVANCED_DETREND - Sáv-specifikus detrending
%
% Inputs:
%   signal - input signal vector
%   method - 'linear', 'polynomial', 'emd'
%   order  - polynomial order (only for polynomial method)
%
% Output:
%   signal_detrended - detrended signal

% Input validation
if any(~isfinite(signal))
    signal = fillmissing(signal, 'nearest');
end

switch lower(method)
    case 'linear'
        % Standard linear detrending
        signal_detrended = detrend(signal, 'linear');
        
    case 'polynomial'
        % Polynomial detrending
        if nargin < 3 || isempty(order)
            order = 2;  % Default quadratic
        end
        
        n = length(signal);
        t = (1:n)';
        
        % Fit polynomial
        p = polyfit(t, signal, order);
        trend = polyval(p, t);
        
        % Remove trend
        signal_detrended = signal - trend;
        
    case 'emd'
        % Empirical Mode Decomposition detrending
        try
            % Check if EMD toolbox is available
            if exist('emd', 'file') == 2
                % Use MATLAB's emd function (R2018a+)
                [imfs, residual] = emd(signal);
                
                % Remove the lowest frequency components (trend)
                % Keep higher frequency IMFs, remove residual
                if size(imfs, 2) > 1
                    signal_detrended = sum(imfs(:, 1:end-1), 2);  % All IMFs except last
                else
                    signal_detrended = signal - residual;  % Remove residual only
                end
                
            else
                % Fallback: use simple EMD implementation
                signal_detrended = simple_emd_detrend(signal);
            end
            
        catch
            % If EMD fails, fallback to polynomial
            warning('EMD failed, using polynomial detrending instead');
            signal_detrended = advanced_detrend(signal, 'polynomial', 3);
        end
        
    otherwise
        % Default to linear
        warning('Unknown detrend method: %s. Using linear.', method);
        signal_detrended = detrend(signal, 'linear');
end

% Ensure output is column vector
signal_detrended = signal_detrended(:);

end

function signal_detrended = simple_emd_detrend(signal)
% Simple EMD implementation for trend removal
% This is a simplified version - for production use proper EMD toolbox

n = length(signal);
t = (1:n)';

% Extract envelope-based trend
upper_env = envelope(signal, 'peak');
lower_env = envelope(signal, 'valley'); 
mean_env = (upper_env + lower_env) / 2;

% Remove trend
signal_detrended = signal - mean_env;

% If envelope function not available, use moving average as fallback
if any(~isfinite(signal_detrended))
    window_size = max(10, round(n/10));
    trend = movmean(signal, window_size);
    signal_detrended = signal - trend;
end

end
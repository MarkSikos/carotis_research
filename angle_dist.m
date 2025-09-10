clear; clc;

%% FEJLETT SZÖGELOSZLÁSOS IDŐSOR ELEMZÉS - IRÁNYMENTI DERIVÁLTFÜGGVÉNNYEL
% Adatelőfeldolgozás + folytonos szögeloszlás elemzés kombinálva

fprintf('=== FEJLETT SZÖGELOSZLÁSOS IDŐSOR ELEMZÉS ===\n');

%% 1. PARAMÉTEREK ÉS BEÁLLÍTÁSOK
data_file = 'data/df_unnormalized.csv';
target_oppart = 4;  % Clamp fázis
min_valid_ratio = 0.70;  % Min 70% valid adat kell

% Szögeloszlás paraméterek
target_fs = 0.5;           % Target mintavételi frekvencia
n_points = 300;            % Interpolált pontok száma (több a finomabb elemzéshez)
target_duration = 600;     % Target időtartam másodpercben

% Deriváltfüggvény paraméterek
derivative_methods = {'gradient', 'local_window', 'savitzky_golay', 'wavelet_based'};
window_sizes = [3, 5, 10, 20];  % Különböző ablakméret a multi-scale elemzéshez
kde_bandwidth = 5;              % Kernel density estimation sávszélesség

% Fiziológiai határok
physio_bounds = [0, 200; 0, 100; 30, 200]; % [SE_min, SE_max; rSO2_min, rSO2_max; MAP_min, MAP_max]

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
    error('Nem találhatóoppart vagy oper_phase oszlop!');
end

patients = unique(phase_data.Identifier);
n_patients = length(patients);
fprintf('Eredeti páciensek száma: %d\n', n_patients);

%% 3. ADATELŐFELDOLGOZÁS ÉS FEJLETT SZÖGELOSZLÁS SZÁMÍTÁS
processed_results = [];
angle_results = struct();
angle_results.patient_data = [];

max_gap_points = round(30 * target_fs);

fprintf('\n=== PÁCIENSENKÉNTI FELDOLGOZÁS ===\n');

for p = 1:n_patients
    patient_id = patients(p);
    patient_data = phase_data(phase_data.Identifier == patient_id, :);
    
    fprintf('Páciens %d feldolgozása...\n', patient_id);
    
    %% 3.1 Adatminőség ellenőrzés (ugyanaz mint előtte)
    SE_signal = patient_data.SE;
    rSO2_signal = patient_data.oper_side_oxig;
    MAP_signal = patient_data.MAP;
    rSO2_contra_signal = patient_data.other_side_oxig;
    
    se_valid_ratio = sum(~isnan(SE_signal)) / length(SE_signal);
    rso2_valid_ratio = sum(~isnan(rSO2_signal)) / length(rSO2_signal);
    map_valid_ratio = sum(~isnan(MAP_signal)) / length(MAP_signal);
    rso2_contra_valid_ratio = sum(~isnan(rSO2_contra_signal)) / length(rSO2_contra_signal);
    
    if se_valid_ratio < min_valid_ratio || rso2_valid_ratio < min_valid_ratio || ...
       map_valid_ratio < min_valid_ratio || rso2_contra_valid_ratio < min_valid_ratio
        fprintf('  ❌ Elégtelen adatminőség\n');
        continue;
    end
    
    try
        %% 3.2-3.5 Előfeldolgozás (mint előtte)
        SE_cleaned = smart_gap_filling(SE_signal, max_gap_points);
        rSO2_cleaned = smart_gap_filling(rSO2_signal, max_gap_points);
        MAP_cleaned = smart_gap_filling(MAP_signal, max_gap_points);
        rSO2_contra_cleaned = smart_gap_filling(rSO2_contra_signal, max_gap_points);
        
        patient_data_for_interp = patient_data;
        patient_data_for_interp.SE = SE_cleaned;
        patient_data_for_interp.oper_side_oxig = rSO2_cleaned;
        patient_data_for_interp.MAP = MAP_cleaned;
        patient_data_for_interp.other_side_oxig = rSO2_contra_cleaned;
        
        interpolated_data = interpolate_data_simple(patient_data_for_interp, n_points, target_duration);
        
        SE_clean = fillmissing(interpolated_data.SE, 'nearest');
        rSO2_clean = fillmissing(interpolated_data.rSO2, 'nearest');
        MAP_clean = fillmissing(interpolated_data.MAP, 'nearest');
        rSO2_contra_clean = fillmissing(interpolated_data.rSO2_contra, 'nearest');
        
        % Fiziológiai ellenőrzés
        SE_physio_valid = SE_clean >= physio_bounds(1,1) & SE_clean <= physio_bounds(1,2);
        rSO2_physio_valid = rSO2_clean >= physio_bounds(2,1) & rSO2_clean <= physio_bounds(2,2);
        MAP_physio_valid = MAP_clean >= physio_bounds(3,1) & MAP_clean <= physio_bounds(3,2);
        rSO2_contra_physio_valid = rSO2_contra_clean >= physio_bounds(2,1) & rSO2_contra_clean <= physio_bounds(2,2);
        
        physio_valid_ratio = sum(SE_physio_valid & rSO2_physio_valid & MAP_physio_valid & rSO2_contra_physio_valid) / length(SE_clean);
        
        if physio_valid_ratio < min_valid_ratio
            fprintf('  ❌ Fiziológiai határok sértése\n');
            continue;
        end
        
        % Detrending és szűrés
        SE_clean = detrend(SE_clean);
        rSO2_clean = detrend(rSO2_clean);
        
        if length(SE_clean) > 20
            cutoff_freq = 0.01;
            [b, a] = butter(1, cutoff_freq/(target_fs/2), 'high');
            SE_clean = filtfilt(b, a, SE_clean);
            rSO2_clean = filtfilt(b, a, rSO2_clean);
        end
        
        %% 4. 🎯 FEJLETT IRÁNYMENTI DERIVÁLTFÜGGVÉNY ELEMZÉS
        
        dt = 1/target_fs; % időlépés
        time_vector = (0:length(SE_clean)-1) * dt;
        
        % 🔥 MULTI-SCALE DERIVÁLT SZÁMÍTÁS
        SE_derivatives = calculate_multiscale_derivatives(SE_clean, dt, window_sizes);
        rSO2_derivatives = calculate_multiscale_derivatives(rSO2_clean, dt, window_sizes);
        
        % 🔥 KONTINUUS SZÖGELOSZLÁS MINDEN SKÁLÁN
        SE_angle_analysis = struct();
        rSO2_angle_analysis = struct();
        
        for scale_idx = 1:length(window_sizes)
            scale_name = sprintf('scale_%d', window_sizes(scale_idx));
            
            % SE szögek ezen a skálán
            SE_angles = atan2d(SE_derivatives.(['scale_' num2str(window_sizes(scale_idx))]), dt);
            rSO2_angles = atan2d(rSO2_derivatives.(['scale_' num2str(window_sizes(scale_idx))]), dt);
            
            % 🔥 KERNEL DENSITY ESTIMATION - folytonos eloszlás
            [SE_kde, SE_kde_x] = calculate_kde(SE_angles, kde_bandwidth);
            [rSO2_kde, rSO2_kde_x] = calculate_kde(rSO2_angles, kde_bandwidth);
            
            % 🔥 FEJLETT STATISZTIKAI ELEMZÉS
            SE_analysis = analyze_continuous_angle_distribution(SE_angles, SE_kde, SE_kde_x);
            rSO2_analysis = analyze_continuous_angle_distribution(rSO2_angles, rSO2_kde, rSO2_kde_x);
            
            % Tárolás
            SE_angle_analysis.(scale_name) = SE_analysis;
            rSO2_angle_analysis.(scale_name) = rSO2_analysis;
        end
        
        % 🔥 WAVELET-ALAPÚ FREKVENCIA-FÜGGŐ SZÖGELOSZLÁS
        if exist('cwt', 'file') == 2
            SE_wavelet_angles = calculate_wavelet_angle_distribution(SE_clean, target_fs);
            rSO2_wavelet_angles = calculate_wavelet_angle_distribution(rSO2_clean, target_fs);
        else
            SE_wavelet_angles = [];
            rSO2_wavelet_angles = [];
        end
        
        % 🔥 KERESZTKORRELÁCIÓ SZÖGELOSZLÁSOK KÖZÖTT
        cross_angle_analysis = analyze_cross_angle_correlations(SE_angle_analysis, rSO2_angle_analysis);
        
        %% 5. EREDMÉNYEK TÁROLÁSA
        result = struct();
        result.PatientID = patient_id;
        result.SE_signal = SE_clean;
        result.rSO2_signal = rSO2_clean;
        result.time_vector = time_vector;
        result.SE_derivatives = SE_derivatives;
        result.rSO2_derivatives = rSO2_derivatives;
        result.SE_angle_analysis = SE_angle_analysis;
        result.rSO2_angle_analysis = rSO2_angle_analysis;
        result.SE_wavelet_angles = SE_wavelet_angles;
        result.rSO2_wavelet_angles = rSO2_wavelet_angles;
        result.cross_angle_analysis = cross_angle_analysis;
        result.data_quality = physio_valid_ratio;
        result.signal_length = length(SE_clean);
        
        angle_results.patient_data = [angle_results.patient_data; result];
        
        fprintf('  ✅ Sikeres (%.1f%% valid, %d multi-scale elemzés)\n', ...
                physio_valid_ratio*100, length(window_sizes));
        
    catch ME
        fprintf('  ❌ Feldolgozási hiba: %s\n', ME.message);
        continue;
    end
end

%% 6. ÖSSZEGZÉS ÉS FEJLETT VIZUALIZÁCIÓ
n_successful = length(angle_results.patient_data);
fprintf('\n=== FELDOLGOZÁS ÖSSZEGZÉS ===\n');
fprintf('Sikeresen feldolgozott páciensek: %d/%d (%.1f%%)\n', ...
        n_successful, n_patients, n_successful/n_patients*100);

if n_successful > 0
    fprintf('\n=== FEJLETT VIZUALIZÁCIÓ ===\n');
    example_patient = angle_results.patient_data(1);
    visualize_advanced_angle_analysis(example_patient, window_sizes);
    
    fprintf('\n=== POPULÁCIÓS MULTI-SCALE ELEMZÉS ===\n');
    analyze_population_multiscale_angles(angle_results.patient_data, window_sizes);
else
    fprintf('❌ Egyetlen páciens sem lett sikeresen feldolgozva!\n');
end

fprintf('\nFejlett szögeloszlásos elemzés befejezve!\n');

%% ============================================================================
%% 🔥 ÚJ FEJLETT HELPER FÜGGVÉNYEK
%% ============================================================================

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

%% KORÁBBI HELPER FÜGGVÉNYEK (változatlanul)
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

function interpolated_data = interpolate_data_simple(patient_data, n_points, target_duration)
    original_time = patient_data.Time;
    original_SE = patient_data.SE;
    original_rSO2 = patient_data.oper_side_oxig;
    original_MAP = patient_data.MAP;
    original_rSO2_contra = patient_data.other_side_oxig;
    
    total_duration = max(original_time) - min(original_time);
    
    if total_duration >= target_duration
        start_time = min(original_time) + (total_duration - target_duration) / 2;
        end_time = start_time + target_duration;
        time_mask = original_time >= start_time & original_time <= end_time;
        selected_time = original_time(time_mask);
        selected_SE = original_SE(time_mask);
        selected_rSO2 = original_rSO2(time_mask);
        selected_MAP = original_MAP(time_mask);
        selected_rSO2_contra = original_rSO2_contra(time_mask);
    else
        selected_time = original_time;
        selected_SE = original_SE;
        selected_rSO2 = original_rSO2;
        selected_MAP = original_MAP;
        selected_rSO2_contra = original_rSO2_contra;
    end
    
    new_time = linspace(min(selected_time), max(selected_time), n_points);
    
    interpolated_data = struct();
    interpolated_data.SE = interp1(selected_time, selected_SE, new_time, 'pchip', 'extrap');
    interpolated_data.rSO2 = interp1(selected_time, selected_rSO2, new_time, 'pchip', 'extrap');
    interpolated_data.MAP = interp1(selected_time, selected_MAP, new_time, 'pchip', 'extrap');
    interpolated_data.rSO2_contra = interp1(selected_time, selected_rSO2_contra, new_time, 'pchip', 'extrap');
    interpolated_data.time = new_time;
end

%% 🎨 FEJLETT VIZUALIZÁCIÓS FÜGGVÉNYEK
function visualize_advanced_angle_analysis(patient_result, window_sizes)
    % Fejlett multi-scale szögeloszlás vizualizáció
    figure('Position', [50, 50, 1600, 1200]);
    
    patient_id = patient_result.PatientID;
    
    % 1. Eredeti jelek és deriváltak
    subplot(4,3,1);
    plot(patient_result.time_vector, patient_result.SE_signal, 'b-', 'LineWidth', 2);
    title(sprintf('SE jel - Páciens %d', patient_id));
    xlabel('Idő (s)'); ylabel('SE');
    grid on;
    
    subplot(4,3,2);
    plot(patient_result.time_vector, patient_result.rSO2_signal, 'r-', 'LineWidth', 2);
    title('rSO2 jel');
    xlabel('Idő (s)'); ylabel('rSO2 (%)');
    grid on;
    
    % 2. Multi-scale deriváltak
    subplot(4,3,3);
    colors = lines(length(window_sizes));
    hold on;
    for i = 1:length(window_sizes)
        scale_name = sprintf('scale_%d', window_sizes(i));
        if isfield(patient_result.SE_derivatives, scale_name)
            deriv = patient_result.SE_derivatives.(scale_name);
            plot(patient_result.time_vector, deriv, 'Color', colors(i,:), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Skála %d', window_sizes(i)));
        end
    end
    title('SE Multi-scale deriváltak');
    xlabel('Idő (s)'); ylabel('dSE/dt');
    legend('Location', 'best');
    grid on;
    
    % 3-6. KDE szögeloszlások különböző skálákon
    for scale_idx = 1:min(4, length(window_sizes))
        subplot(4,3,3+scale_idx);
        scale_name = sprintf('scale_%d', window_sizes(scale_idx));
        
        if isfield(patient_result.SE_angle_analysis, scale_name)
            SE_analysis = patient_result.SE_angle_analysis.(scale_name);
            rSO2_analysis = patient_result.rSO2_angle_analysis.(scale_name);
            
            % KDE plotok
            plot(SE_analysis.kde_x, SE_analysis.kde_values, 'b-', 'LineWidth', 2, 'DisplayName', 'SE');
            hold on;
            plot(rSO2_analysis.kde_x, rSO2_analysis.kde_values, 'r-', 'LineWidth', 2, 'DisplayName', 'rSO2');
            
            % Csúcsok jelölése
            if SE_analysis.n_peaks > 0
                plot(SE_analysis.peak_angles, SE_analysis.peak_heights, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            end
            if rSO2_analysis.n_peaks > 0
                plot(rSO2_analysis.peak_angles, rSO2_analysis.peak_heights, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            end
            
            xline(0, '--k');
            title(sprintf('Szögeloszlás - Skála %d', window_sizes(scale_idx)));
            xlabel('Szög (°)'); ylabel('Sűrűség');
            legend('Location', 'best');
            grid on;
        end
    end
    
    % 7. Keresztkorreláció heatmap
    subplot(4,3,8);
    correlations = [];
    scale_labels = {};
    for i = 1:length(window_sizes)
        scale_name = sprintf('scale_%d', window_sizes(i));
        if isfield(patient_result.cross_angle_analysis, scale_name)
            correlations(end+1) = patient_result.cross_angle_analysis.(scale_name).correlation;
            scale_labels{end+1} = sprintf('S%d', window_sizes(i));
        end
    end
    
    if ~isempty(correlations)
        bar(1:length(correlations), correlations, 'FaceColor', [0.7 0.3 0.9]);
        set(gca, 'XTick', 1:length(correlations), 'XTickLabel', scale_labels);
        title('SE-rSO2 Keresztkorreláció');
        ylabel('Korreláció');
        grid on;
        ylim([-1 1]);
    end
    
    % 8. Frekvencia-függő wavelet elemzés
    subplot(4,3,9);
    if ~isempty(patient_result.SE_wavelet_angles)
        band_names = fieldnames(patient_result.SE_wavelet_angles);
        for i = 1:length(band_names)
            band_angles = patient_result.SE_wavelet_angles.(band_names{i});
            if ~isempty(band_angles)
                histogram(band_angles, 20, 'Normalization', 'pdf', ...
                         'DisplayName', band_names{i});
                hold on;
            end
        end
        title('Wavelet frekvenciasáv szögek');
        xlabel('Szög (°)'); ylabel('Sűrűség');
        legend('Location', 'best');
        grid on;
    else
        text(0.5, 0.5, 'Wavelet elemzés nem elérhető', 'HorizontalAlignment', 'center');
        title('Wavelet elemzés');
    end
    
    % 9. Irányspecifikus statisztikák
    subplot(4,3,10);
    if length(window_sizes) >= 1
        scale_name = sprintf('scale_%d', window_sizes(1));
        if isfield(patient_result.SE_angle_analysis, scale_name)
            SE_stats = patient_result.SE_angle_analysis.(scale_name);
            rSO2_stats = patient_result.rSO2_angle_analysis.(scale_name);
            
            categories = {'Felfelé%', 'Lefelé%', 'Entropy', 'Aszimmetria'};
            SE_vals = [SE_stats.positive_ratio*100, SE_stats.negative_ratio*100, ...
                      SE_stats.kde_entropy, SE_stats.directional_asymmetry];
            rSO2_vals = [rSO2_stats.positive_ratio*100, rSO2_stats.negative_ratio*100, ...
                        rSO2_stats.kde_entropy, rSO2_stats.directional_asymmetry];
            
            x = 1:length(categories);
            bar(x-0.2, SE_vals, 0.4, 'FaceColor', [0.2 0.6 0.8], 'DisplayName', 'SE');
            hold on;
            bar(x+0.2, rSO2_vals, 0.4, 'FaceColor', [0.8 0.4 0.2], 'DisplayName', 'rSO2');
            
            set(gca, 'XTick', x, 'XTickLabel', categories);
            title('Irányspecifikus metrikák');
            legend('Location', 'best');
            grid on;
        end
    end
    
    % 10-12. További analitikai plotok...
    subplot(4,3,11);
    text(0.1, 0.8, sprintf('Páciens ID: %d', patient_id), 'FontSize', 12);
    text(0.1, 0.7, sprintf('Adatminőség: %.1f%%', patient_result.data_quality*100), 'FontSize', 12);
    text(0.1, 0.6, sprintf('Jelhossz: %d pont', patient_result.signal_length), 'FontSize', 12);
    text(0.1, 0.5, sprintf('Elemzett skálák: %d', length(window_sizes)), 'FontSize', 12);
    axis off;
    title('Összefoglaló');
    
    sgtitle(sprintf('Fejlett Multi-scale Szögeloszlás Elemzés - Páciens %d', patient_id), ...
            'FontSize', 16, 'FontWeight', 'bold');
end

function analyze_population_multiscale_angles(patient_data_array, window_sizes)
    % Populációs szintű multi-scale elemzés
    n_patients = length(patient_data_array);
    fprintf('Multi-scale populációs elemzés %d páciensen:\n', n_patients);
    
    for scale_idx = 1:length(window_sizes)
        scale_name = sprintf('scale_%d', window_sizes(scale_idx));
        fprintf('\n--- Skála %d ---\n', window_sizes(scale_idx));
        
        all_SE_entropy = [];
        all_rSO2_entropy = [];
        all_correlations = [];
        all_SE_asymmetry = [];
        all_rSO2_asymmetry = [];
        
        for i = 1:n_patients
            patient = patient_data_array(i);
            
            if isfield(patient.SE_angle_analysis, scale_name)
                SE_stats = patient.SE_angle_analysis.(scale_name);
                rSO2_stats = patient.rSO2_angle_analysis.(scale_name);
                
                all_SE_entropy(end+1) = SE_stats.kde_entropy;
                all_rSO2_entropy(end+1) = rSO2_stats.kde_entropy;
                all_SE_asymmetry(end+1) = SE_stats.directional_asymmetry;
                all_rSO2_asymmetry(end+1) = rSO2_stats.directional_asymmetry;
                
                if isfield(patient.cross_angle_analysis, scale_name)
                    all_correlations(end+1) = patient.cross_angle_analysis.(scale_name).correlation;
                end
            end
        end
        
        if ~isempty(all_SE_entropy)
            fprintf('SE entropy: %.3f ± %.3f\n', mean(all_SE_entropy), std(all_SE_entropy));
            fprintf('rSO2 entropy: %.3f ± %.3f\n', mean(all_rSO2_entropy), std(all_rSO2_entropy));
            fprintf('SE aszimmetria: %.3f ± %.3f\n', mean(all_SE_asymmetry), std(all_SE_asymmetry));
            fprintf('rSO2 aszimmetria: %.3f ± %.3f\n', mean(all_rSO2_asymmetry), std(all_rSO2_asymmetry));
            
            if ~isempty(all_correlations)
                fprintf('Átlagos SE-rSO2 korreláció: %.3f ± %.3f\n', mean(all_correlations), std(all_correlations));
            end
        end
    end
end
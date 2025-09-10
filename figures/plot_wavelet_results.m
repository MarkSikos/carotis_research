function plot_wavelet_results(csv_file)
% PLOT_BAND_SPECIFIC_RESULTS - Sáv-specifikus eredmények vizualizálása
%
% Használat:
%   plot_band_specific_results('band_specific_results.csv');

%% 1. Adatok betöltése
fprintf('=== BAND-SPECIFIC RESULTS PLOTTING ===\n');

if ~exist(csv_file, 'file')
    error('Fájl nem található: %s', csv_file);
end

% CSV beolvasása
data = readtable(csv_file);
fprintf('Adatok betöltve: %s\n', csv_file);

% Oszlopnevek ellenőrzése
fprintf('Elérhető oszlopok: ');
disp(data.Properties.VariableNames);

%% 2. Adatok kinyerése
band_names = data.Band;
if iscell(band_names)
    % Cell array stringekké alakítása
    band_names = cellfun(@(x) x, band_names, 'UniformOutput', false);
else
    band_names = cellstr(band_names);
end

coherence_values = data.Coherence;
PLV_values = data.PLV;
phase_deg = data.CMP_deg;
COx_values = data.COx;
n_patients = data.N_patients;

n_bands = length(band_names);
fprintf('Betöltött sávok: %d\n', n_bands);

%% 3. Színek definiálása
colors = [0.8, 0.2, 0.2;    % Endothelial - Piros
          0.2, 0.8, 0.2;    % Neurogenic - Zöld
          0.2, 0.2, 0.8;    % Myogenic - Kék
          0.8, 0.8, 0.2;    % Respiratory - Sárga
          0.8, 0.2, 0.8];   % Cardiac - Magenta

% Ha több sáv van, mint szín, generáljunk többet
if n_bands > size(colors, 1)
    extra_colors = hsv(n_bands - size(colors, 1));
    colors = [colors; extra_colors];
end

%% 4. PLOT 1: MAIN COMPARISON DASHBOARD
figure('Position', [100, 100, 1600, 900], 'Name', 'Band-Specific Results Dashboard');

% Subplot 1: Coherence comparison
subplot(2, 4, 1);
for i = 1:n_bands
    bar(i, coherence_values(i), 'FaceColor', colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean Coherence');
title('Coherence by Frequency Band');
xtickangle(45);
grid on;
ylim([0, max(coherence_values)*1.2]);

% Értékek kiírása
for i = 1:n_bands
    text(i, coherence_values(i) + max(coherence_values)*0.05, ...
         sprintf('%.3f', coherence_values(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

% Subplot 2: PLV comparison
subplot(2, 4, 2);
for i = 1:n_bands
    bar(i, PLV_values(i), 'FaceColor', colors(i,:)*0.7, 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean PLV');
title('Phase Locking Value by Band');
xtickangle(45);
grid on;
ylim([0, max(PLV_values)*1.2]);

% Értékek kiírása
for i = 1:n_bands
    text(i, PLV_values(i) + max(PLV_values)*0.05, ...
         sprintf('%.3f', PLV_values(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

% Subplot 3: Phase Direction
subplot(2, 4, 3);
for i = 1:n_bands
    color_intensity = 0.8 + 0.2 * (phase_deg(i) < 0);  % Negatív fázisok sötétebbek
    bar(i, phase_deg(i), 'FaceColor', colors(i,:)*color_intensity, ...
        'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean Phase (degrees)');
title('Phase Direction by Band');
xtickangle(45);
grid on;
line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
ylim([min(phase_deg)*1.2, max(phase_deg)*1.2]);

% Értékek kiírása
for i = 1:n_bands
    y_offset = sign(phase_deg(i)) * abs(max(phase_deg)-min(phase_deg)) * 0.05;
    text(i, phase_deg(i) + y_offset, sprintf('%.0f°', phase_deg(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

% Subplot 4: COx Validation
subplot(2, 4, 4);
for i = 1:n_bands
    if ~isnan(COx_values(i))
        color_intensity = 0.5 + 0.5 * (COx_values(i) > 0);  % Pozitív értékek világosabbak
        bar(i, COx_values(i), 'FaceColor', colors(i,:)*color_intensity, ...
            'EdgeColor', 'k', 'LineWidth', 1.5);
        hold on;
    end
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean COx');
title('COx Validation by Band');
xtickangle(45);
grid on;
line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);

% Értékek kiírása
for i = 1:n_bands
    if ~isnan(COx_values(i))
        y_offset = sign(COx_values(i)) * abs(max(COx_values)-min(COx_values)) * 0.1;
        text(i, COx_values(i) + y_offset, sprintf('%.3f', COx_values(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
    end
end
% Subplot 5-6: Combined Coherence + PLV (kompaktabb)
subplot(2, 4, 5);
x_pos = 1:n_bands;
bar_width = 0.35;

% Coherence bars
for i = 1:n_bands
    bar(x_pos(i) - bar_width/2, coherence_values(i), bar_width, ...
        'FaceColor', colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1);
    hold on;
end

% PLV bars (scaled to same range)
plv_scaled = PLV_values * max(coherence_values) / max(PLV_values);
for i = 1:n_bands
    bar(x_pos(i) + bar_width/2, plv_scaled(i), bar_width, ...
        'FaceColor', colors(i,:)*0.6, 'EdgeColor', 'k', 'LineWidth', 1);
end

set(gca, 'XTick', x_pos, 'XTickLabel', band_names);
xlabel('Frequency Bands');
ylabel('Values');
title('Coherence vs PLV');
xtickangle(45);
grid on;
legend({'Coherence', 'PLV'}, 'Location', 'best', 'FontSize', 8);

% Subplot 6: Polar Phase Plot (kompaktabb)
subplot(2, 4, 6);
phase_rad = deg2rad(phase_deg);

polarplot(phase_rad, abs(phase_rad), 'o', 'MarkerSize', 8, 'LineWidth', 2);
hold on;

for i = 1:n_bands
    polarplot([phase_rad(i), phase_rad(i)], [0, abs(phase_rad(i))], ...
              'Color', colors(i,:), 'LineWidth', 3);
end

% Reference lines
polarplot([0, 2*pi], [pi/2, pi/2], 'k--');
polarplot([0, 2*pi], [pi, pi], 'k--');

title('Phase Directions', 'FontSize', 10);

% Subplot 7: Quick Stats Table
subplot(2, 4, 7);
axis off;

% Compact summary table
stats_text = {};
stats_text{1} = 'SUMMARY STATS';
stats_text{2} = '=============';
stats_text{3} = '';

% Find best performers
[max_coh, max_coh_idx] = max(coherence_values);
[max_plv, max_plv_idx] = max(PLV_values);

stats_text{4} = sprintf('Best Coherence:');
stats_text{5} = sprintf('  %s (%.3f)', band_names{max_coh_idx}, max_coh);
stats_text{6} = '';
stats_text{7} = sprintf('Best PLV:');
stats_text{8} = sprintf('  %s (%.3f)', band_names{max_plv_idx}, max_plv);
stats_text{9} = '';
stats_text{10} = sprintf('Overall:');
stats_text{11} = sprintf('  Mean Coh: %.3f', mean(coherence_values));
stats_text{12} = sprintf('  Mean PLV: %.3f', mean(PLV_values));

overall_coh = mean(coherence_values);
overall_plv = mean(PLV_values);

if overall_coh > 0.4 && overall_plv > 0.6
    assessment = 'EXCELLENT';
    assess_color = [0.2, 0.8, 0.2];
elseif overall_coh > 0.3 && overall_plv > 0.5
    assessment = 'GOOD';
    assess_color = [0.8, 0.8, 0.2];
elseif overall_coh > 0.2 && overall_plv > 0.3
    assessment = 'MODERATE';
    assess_color = [0.8, 0.6, 0.2];
else
    assessment = 'WEAK';
    assess_color = [0.8, 0.2, 0.2];
end

text(0.05, 0.95, stats_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontName', 'FixedWidth', ...
     'FontSize', 12, 'BackgroundColor', assess_color, ...
     'EdgeColor', [0.3, 0.3, 0.3], 'LineWidth', 2);
% PLOT_BAND_SPECIFIC_RESULTS - Sáv-specifikus eredmények vizualizálása
%
% Használat:
%   plot_band_specific_results('band_specific_results.csv');

%% 1. Adatok betöltése
fprintf('=== BAND-SPECIFIC RESULTS PLOTTING ===\n');

if ~exist(csv_file, 'file')
    error('Fájl nem található: %s', csv_file);
end

% CSV beolvasása
data = readtable(csv_file);
fprintf('Adatok betöltve: %s\n', csv_file);

% Oszlopnevek ellenőrzése
fprintf('Elérhető oszlopok: ');
disp(data.Properties.VariableNames);

%% 2. Adatok kinyerése
band_names = data.Band;
if iscell(band_names)
    % Cell array stringekké alakítása
    band_names = cellfun(@(x) x, band_names, 'UniformOutput', false);
else
    band_names = cellstr(band_names);
end

coherence_values = data.Coherence;
PLV_values = data.PLV;
phase_deg = data.CMP_deg;
COx_values = data.COx;
n_patients = data.N_patients;

n_bands = length(band_names);
fprintf('Betöltött sávok: %d\n', n_bands);

%% 3. Színek definiálása
colors = [0.8, 0.2, 0.2;    % Endothelial - Piros
          0.2, 0.8, 0.2;    % Neurogenic - Zöld
          0.2, 0.2, 0.8;    % Myogenic - Kék
          0.8, 0.8, 0.2;    % Respiratory - Sárga
          0.8, 0.2, 0.8];   % Cardiac - Magenta

% Ha több sáv van, mint szín, generáljunk többet
if n_bands > size(colors, 1)
    extra_colors = hsv(n_bands - size(colors, 1));
    colors = [colors; extra_colors];
end

%% 4. PLOT 1: COMPACT COMPARISON DASHBOARD  
figure('Position', [100, 100, 1400, 700], 'Name', 'Band-Specific Results Dashboard');

% Subplot 1: Coherence comparison
subplot(2, 4, 1);
for i = 1:n_bands
    bar(i, coherence_values(i), 'FaceColor', colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean Coherence');
title('Coherence by Frequency Band');
xtickangle(45);
grid on;
ylim([0, max(coherence_values)*1.2]);

% Értékek kiírása
for i = 1:n_bands
    text(i, coherence_values(i) + max(coherence_values)*0.05, ...
         sprintf('%.3f', coherence_values(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

% Subplot 2: PLV comparison
subplot(2, 4, 2);
for i = 1:n_bands
    bar(i, PLV_values(i), 'FaceColor', colors(i,:)*0.7, 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean PLV');
title('Phase Locking Value by Band');
xtickangle(45);
grid on;
ylim([0, max(PLV_values)*1.2]);

% Értékek kiírása
for i = 1:n_bands
    text(i, PLV_values(i) + max(PLV_values)*0.05, ...
         sprintf('%.3f', PLV_values(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

% Subplot 3: Phase Direction
subplot(2, 4, 3);
for i = 1:n_bands
    color_intensity = 0.8 + 0.2 * (phase_deg(i) < 0);  % Negatív fázisok sötétebbek
    bar(i, phase_deg(i), 'FaceColor', colors(i,:)*color_intensity, ...
        'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean Phase (degrees)');
title('Phase Direction by Band');
xtickangle(45);
grid on;
line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
ylim([min(phase_deg)*1.2, max(phase_deg)*1.2]);

% Értékek kiírása
for i = 1:n_bands
    y_offset = sign(phase_deg(i)) * abs(max(phase_deg)-min(phase_deg)) * 0.05;
    text(i, phase_deg(i) + y_offset, sprintf('%.0f°', phase_deg(i)), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

% Subplot 4: COx Validation
subplot(2, 4, 4);
for i = 1:n_bands
    if ~isnan(COx_values(i))
        color_intensity = 0.5 + 0.5 * (COx_values(i) > 0);  % Pozitív értékek világosabbak
        bar(i, COx_values(i), 'FaceColor', colors(i,:)*color_intensity, ...
            'EdgeColor', 'k', 'LineWidth', 1.5);
        hold on;
    end
end
set(gca, 'XTick', 1:n_bands, 'XTickLabel', band_names);
ylabel('Mean COx');
title('COx Validation by Band');
xtickangle(45);
grid on;
line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);

% Értékek kiírása
for i = 1:n_bands
    if ~isnan(COx_values(i))
        y_offset = sign(COx_values(i)) * abs(max(COx_values)-min(COx_values)) * 0.1;
        text(i, COx_values(i) + y_offset, sprintf('%.3f', COx_values(i)), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
    end
end

% Subplot 5-6: Combined Coherence + PLV
subplot(2, 4, [5, 6]);
x_pos = 1:n_bands;
bar_width = 0.35;

% Coherence bars
for i = 1:n_bands
    bar(x_pos(i) - bar_width/2, coherence_values(i), bar_width, ...
        'FaceColor', colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1);
    hold on;
end

% PLV bars (scaled to same range)
plv_scaled = PLV_values * max(coherence_values) / max(PLV_values);
for i = 1:n_bands
    bar(x_pos(i) + bar_width/2, plv_scaled(i), bar_width, ...
        'FaceColor', colors(i,:)*0.6, 'EdgeColor', 'k', 'LineWidth', 1);
end

set(gca, 'XTick', x_pos, 'XTickLabel', band_names);
xlabel('Frequency Bands');
ylabel('Normalized Values');
title('Coherence vs PLV Comparison');
xtickangle(45);
grid on;

% Dual y-axis labels
yyaxis left
ylabel('Coherence');
ylim([0, max(coherence_values)*1.2]);

yyaxis right
ylabel('PLV');
ylim([0, max(PLV_values)*1.2]);

legend({'Coherence', 'PLV (scaled)'}, 'Location', 'best');

% Subplot 7: Polar Phase Plot
subplot(2, 4, 7);
phase_rad = deg2rad(phase_deg);
theta = linspace(0, 2*pi, n_bands+1);
theta = theta(1:end-1);  % Remove duplicate point

polarplot(phase_rad, abs(phase_rad), 'o', 'MarkerSize', 10, 'LineWidth', 2);
hold on;

for i = 1:n_bands
    polarplot([phase_rad(i), phase_rad(i)], [0, abs(phase_rad(i))], ...
              'Color', colors(i,:), 'LineWidth', 4);
    
    % Band labels
    [x, y] = pol2cart(phase_rad(i), abs(phase_rad(i))*1.2);
    text(x, y, band_names{i}, 'HorizontalAlignment', 'center', ...
         'FontSize', 8, 'FontWeight', 'bold');
end

% Reference lines
polarplot([0, 2*pi], [pi/2, pi/2], 'k--');
polarplot([0, 2*pi], [pi, pi], 'k--');

title('Phase Directions');

% Subplot 8: Statistics Summary
subplot(2, 4, 8);
axis off;

% Create summary text
summary_text = {};
summary_text{1} = 'BAND STATISTICS SUMMARY';
summary_text{2} = '========================';
summary_text{3} = '';

for i = 1:n_bands
    summary_text{end+1} = sprintf('%s Band:', band_names{i});
    summary_text{end+1} = sprintf('  Coherence: %.3f', coherence_values(i));
    summary_text{end+1} = sprintf('  PLV: %.3f', PLV_values(i));
    summary_text{end+1} = sprintf('  Phase: %.0f°', phase_deg(i));
    summary_text{end+1} = sprintf('  COx: %.3f', COx_values(i));
    summary_text{end+1} = sprintf('  N patients: %d', n_patients(i));
    summary_text{end+1} = '';
end

% Overall assessment
overall_coh = mean(coherence_values);
overall_plv = mean(PLV_values);

summary_text{end+1} = 'OVERALL ASSESSMENT:';
summary_text{end+1} = sprintf('  Mean Coherence: %.3f', overall_coh);
summary_text{end+1} = sprintf('  Mean PLV: %.3f', overall_plv);

if overall_coh > 0.4 && overall_plv > 0.6
    assessment = 'EXCELLENT';
elseif overall_coh > 0.3 && overall_plv > 0.5
    assessment = 'GOOD';
elseif overall_coh > 0.2 && overall_plv > 0.3
    assessment = 'MODERATE';
else
    assessment = 'WEAK';
end

summary_text{end+1} = sprintf('  Coupling: %s', assessment);

text(0.05, 0.95, summary_text, 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontName', 'FixedWidth', ...
     'FontSize', 9, 'BackgroundColor', [0.95, 0.95, 0.95], ...
     'EdgeColor', [0.5, 0.5, 0.5]);

sgtitle('Band-Specific Wavelet Coupling Analysis - Comprehensive Results');

%%
% Console kiírás a subplot helyett
fprintf('\n=== DETAILED BAND STATISTICS ===\n');
for i = 1:n_bands
    fprintf('%s Band:\n', band_names{i});
    fprintf('  Coherence: %.3f\n', coherence_values(i));
    fprintf('  PLV: %.3f\n', PLV_values(i));
    fprintf('  Phase: %.0f°\n', phase_deg(i));
    fprintf('  COx: %.3f\n', COx_values(i));
    fprintf('  N patients: %d\n', n_patients(i));
    fprintf('\n');
end

%% 6. PLOT 3: DETAILED HEATMAP
figure('Position', [200, 200, 1000, 600], 'Name', 'Band Results Heatmap');

% Prepare data for heatmap
metrics = {'Coherence', 'PLV', '|Phase|°', '|COx|', 'N_patients'};
heatmap_data = [coherence_values, PLV_values, abs(phase_deg), abs(COx_values), n_patients/max(n_patients)];

% Create heatmap
h = heatmap(band_names, metrics, heatmap_data', 'Colormap', parula, 'ColorbarVisible', 'on');
h.Title = 'Band-Specific Results Heatmap';
h.XLabel = 'Frequency Bands';
h.YLabel = 'Metrics';

% Add text annotations
for i = 1:length(band_names)
    for j = 1:length(metrics)
        if j == 5  % N_patients - show actual numbers
            value_text = sprintf('%d', n_patients(i));
        elseif j == 3  % Phase - show degrees
            value_text = sprintf('%.0f°', phase_deg(i));
        else
            value_text = sprintf('%.3f', heatmap_data(i, j));
        end
        
        % Add text (Note: this might not work in all MATLAB versions)
        try
            h.NodeChildren(3).NodeChildren(j).NodeChildren(i).Text = value_text;
        catch
            % If text annotation fails, skip
        end
    end
end

%% 7. Console Report
fprintf('\n=== BAND-SPECIFIC PLOTTING COMPLETE ===\n');
fprintf('Processed bands: %d\n', n_bands);
fprintf('Generated plots: Basic Metrics, Comparison & Summary, Heatmap\n');

fprintf('\nRESULTS SUMMARY:\n');
[~, best_coh_idx] = max(coherence_values);
[~, best_plv_idx] = max(PLV_values);

fprintf('Best Coherence: %s (%.3f)\n', band_names{best_coh_idx}, coherence_values(best_coh_idx));
fprintf('Best PLV: %s (%.3f)\n', band_names{best_plv_idx}, PLV_values(best_plv_idx));
fprintf('Overall Assessment: %s coupling\n', assessment);

% Coupling interpretation
fprintf('\nCOUPLING INTERPRETATION:\n');
for i = 1:n_bands
    if phase_deg(i) > 0
        direction = 'SE leads rSO2';
    elseif phase_deg(i) < -90
        direction = 'rSO2 strongly leads SE';
    elseif phase_deg(i) < 0
        direction = 'rSO2 leads SE';
    else
        direction = 'Synchronized';
    end
    
    fprintf('  %s: %s (%.0f°)\n', band_names{i}, direction, phase_deg(i));
end

end
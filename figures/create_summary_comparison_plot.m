
function create_summary_comparison_plot(results_all_bands, band_names)
% Összehasonlító plot
figure('Name', 'Band Comparison', 'Position', [200, 200, 1400, 600]);

valid_bands = {};
coherence_vals = [];
PLV_vals = [];
CMP_vals = [];

for i = 1:length(band_names)
    if isfield(results_all_bands, band_names{i}) && ~isempty(results_all_bands.(band_names{i}))
        valid_bands{end+1} = band_names{i};
        result = results_all_bands.(band_names{i});
        coherence_vals(end+1) = result.mean_coherence;
        PLV_vals(end+1) = result.mean_PLV;
        CMP_vals(end+1) = rad2deg(result.mean_CMP);
    end
end

% Subplot 1: Coherence comparison
subplot(1, 3, 1);
bar(coherence_vals);
set(gca, 'XTickLabel', valid_bands);
ylabel('Mean Coherence');
title('Coherence Comparison');
xtickangle(45);

% Subplot 2: PLV comparison
subplot(1, 3, 2);
bar(PLV_vals);
set(gca, 'XTickLabel', valid_bands);
ylabel('Mean PLV');
title('PLV Comparison');
xtickangle(45);

% Subplot 3: CMP comparison
subplot(1, 3, 3);
bar(CMP_vals);
set(gca, 'XTickLabel', valid_bands);
ylabel('Mean CMP (degrees)');
title('Phase Direction Comparison');
xtickangle(45);
line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');

sgtitle('Frequency Band Comparison');
savefig('band_comparison.fig');
end
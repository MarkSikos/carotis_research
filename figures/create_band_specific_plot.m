
function create_band_specific_plot(band_results, band_name, band_params)
% Sáv-specifikus plot készítése
figure('Name', sprintf('%s Band Analysis', band_name), 'Position', [100, 100, 1200, 800]);

% Subplot 1: Coherence distribution
subplot(2, 3, 1);
histogram(band_results.coherence_values, 'BinWidth', 0.05);
xlabel('Coherence');
ylabel('Páciensek száma');
title(sprintf('%s Coherence Distribution', band_name));
line([band_results.mean_coherence, band_results.mean_coherence], ylim, 'Color', 'r', 'LineWidth', 2);

% Subplot 2: PLV distribution  
subplot(2, 3, 2);
histogram(band_results.PLV_values, 'BinWidth', 0.02);
xlabel('PLV');
ylabel('Páciensek száma');
title(sprintf('%s PLV Distribution', band_name));
line([band_results.mean_PLV, band_results.mean_PLV], ylim, 'Color', 'r', 'LineWidth', 2);

% Subplot 3: CMP polar plot
subplot(2, 3, 3);
polarhistogram(band_results.CMP_values, 'BinWidth', pi/8);
title(sprintf('%s Circular Mean Phase', band_name));

% Subplot 4: COx validation
subplot(2, 3, 4);
if ~all(isnan(band_results.COx_values))
    histogram(band_results.COx_values(~isnan(band_results.COx_values)), 'BinWidth', 0.1);
    xlabel('COx');
    ylabel('Páciensek száma');
    title(sprintf('%s COx Validation', band_name));
    line([band_results.mean_COx, band_results.mean_COx], ylim, 'Color', 'r', 'LineWidth', 2);
end

% Subplot 5: Coherence vs PLV correlation
subplot(2, 3, 5);
scatter(band_results.coherence_values, band_results.PLV_values);
xlabel('Coherence');
ylabel('PLV');
title(sprintf('%s: Coherence vs PLV', band_name));
lsline;

% Subplot 6: Summary stats
subplot(2, 3, 6);
stats = [band_results.mean_coherence; band_results.mean_PLV; abs(band_results.mean_COx)];
bar(stats);
set(gca, 'XTickLabel', {'Coherence', 'PLV', '|COx|'});
ylabel('Érték');
title(sprintf('%s Summary Stats', band_name));

sgtitle(sprintf('%s Band Analysis (%.3f-%.3f Hz, N=%d)', ...
               band_name, band_params.range(1), band_params.range(2), ...
               length(band_results.valid_patients)));

% Save
savefig(sprintf('%s_band_analysis.fig', band_name));
end

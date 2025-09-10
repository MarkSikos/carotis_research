%% CORRELATION ANALYSIS SCRIPT
%% ========================================================================
% Cerebrális autoregulációs metrikák korrelációs analízise
% Input: specific_indices_matrix.csv
% Output: Korrelációs mátrix, heatmap, top korreláció párok
%% ========================================================================

function correlation_analysis_results = analyze_metric_correlations()
    fprintf('=== CEREBRÁLIS AUTOREGULÁCIÓS METRIKÁK KORRELÁCIÓS ANALÍZISE ===\n\n');
    
    % CSV beolvasása
    try
        data_table = readtable('specific_indices_matrix.csv');
        fprintf('CSV sikeresen beolvasva: %d páciensek × %d változók\n', height(data_table), width(data_table));
    catch
        error('Hiba a CSV beolvasásakor. Ellenőrizd hogy a specific_indices_matrix.csv létezik!');
    end
    
    % PatientID oszlop eltávolítása
    if ismember('PatientID', data_table.Properties.VariableNames)
        data_numeric = data_table(:, 2:end); % PatientID kihagyása
        fprintf('PatientID oszlop eltávolítva. Elemzendő metrikák: %d\n', width(data_numeric));
    else
        data_numeric = data_table;
    end
    
    % Numerikus mátrix konvertálás
    data_matrix = table2array(data_numeric);
    metric_names = data_numeric.Properties.VariableNames;
    
    % NaN értékek ellenőrzése
    nan_counts = sum(isnan(data_matrix), 1);
    valid_metrics = nan_counts < height(data_table) * 0.8; % Max 80% NaN
    
    fprintf('Metrikák NaN eloszlása:\n');
    fprintf('  Teljes metrikák: %d\n', sum(valid_metrics));
    fprintf('  Túl sok NaN (>80%%): %d\n', sum(~valid_metrics));
    
    % Csak valid metrikák megtartása
    data_clean = data_matrix(:, valid_metrics);
    metrics_clean = metric_names(valid_metrics);
    
    fprintf('\nFinals analízis: %d metrika\n', length(metrics_clean));
    
    % Korrelációs mátrix számítása
    fprintf('\nKorrelációs mátrix számítása...\n');
    correlation_matrix = corr(data_clean, 'rows', 'pairwise');
    
    % Eredmények struct
    results = struct();
    results.data_table = data_table;
    results.correlation_matrix = correlation_matrix;
    results.metric_names = metrics_clean;
    results.data_clean = data_clean;
    
    % Top korrelációk keresése
    fprintf('Top korrelációk elemzése...\n');
    top_correlations = find_top_correlations(correlation_matrix, metrics_clean);
    results.top_correlations = top_correlations;
    
    % Heatmap generálása
    fprintf('Korrelációs heatmap generálása...\n');
    generate_correlation_heatmap(correlation_matrix, metrics_clean);
    
    % Sáv-specifikus korrelációs analízis
    fprintf('Sáv-specifikus korrelációs analízis...\n');
    band_analysis = analyze_band_correlations(correlation_matrix, metrics_clean);
    results.band_analysis = band_analysis;
    
    % Autoregulációs metrikák korrelációja
    fprintf('Autoregulációs metrikák korrelációja...\n');
    autoreg_analysis = analyze_autoregulation_correlations(correlation_matrix, metrics_clean);
    results.autoregulation_analysis = autoreg_analysis;
    
    % Eredmények mentése
    save('correlation_analysis_results.mat', 'results');
    fprintf('\n=== KORRELÁCIÓS ANALÍZIS BEFEJEZVE ===\n');
    fprintf('Eredmények mentve: correlation_analysis_results.mat\n');
    
    correlation_analysis_results = results;
end

function top_corr = find_top_correlations(corr_matrix, metric_names)
    % Top 50 korreláció keresése (saját magával való korrelációk nélkül)
    n = size(corr_matrix, 1);
    
    % Upper triangle maszk (saját magával való korrelációk nélkül)
    mask = triu(true(n), 1);
    
    % Korrelációs értékek és indexek
    corr_values = corr_matrix(mask);
    [row_idx, col_idx] = find(mask);
    
    % Abszolút értékek szerint rendezés
    [sorted_corr, sort_idx] = sort(abs(corr_values), 'descend');
    
    % Top 50 (vagy kevesebb ha nincs elég)
    top_n = min(50, length(sorted_corr));
    top_indices = sort_idx(1:top_n);
    
    % Eredmények táblázat
    top_corr = table();
    top_corr.Metric1 = metric_names(row_idx(top_indices))';
    top_corr.Metric2 = metric_names(col_idx(top_indices))';
    top_corr.Correlation = corr_values(top_indices);
    top_corr.AbsCorrelation = abs(corr_values(top_indices));
    
    fprintf('\nTOP 10 LEGERŐSEBB KORRELÁCIÓ:\n');
    for i = 1:min(10, height(top_corr))
        fprintf('%d. %s <-> %s: r = %.4f\n', i, ...
            top_corr.Metric1{i}, top_corr.Metric2{i}, top_corr.Correlation(i));
    end
    
    % CSV export
    writetable(top_corr, 'top_correlations.csv');
    fprintf('\nTop korrelációk exportálva: top_correlations.csv\n');
end

function generate_correlation_heatmap(corr_matrix, metric_names)
    % Korrelációs heatmap generálása
    
    % Clustered heatmap (ha sok metrika van)
    if length(metric_names) > 50
        % Hierarchikus klaszterezés
        dist_matrix = 1 - abs(corr_matrix);
        linkage_matrix = linkage(squareform(dist_matrix), 'average');
        cluster_order = optimalleaforder(linkage_matrix, dist_matrix);
        
        % Átrendezett mátrix
        corr_clustered = corr_matrix(cluster_order, cluster_order);
        names_clustered = metric_names(cluster_order);
    else
        corr_clustered = corr_matrix;
        names_clustered = metric_names;
    end
    
    % Heatmap plot
    figure('Position', [100, 100, 1200, 1000]);
    
    % Custom colormap (blue-white-red)
    n_colors = 256;
    colors = [linspace(0, 1, n_colors/2)', linspace(0, 1, n_colors/2)', ones(n_colors/2, 1); ...
              ones(n_colors/2, 1), linspace(1, 0, n_colors/2)', linspace(1, 0, n_colors/2)'];
    
    imagesc(corr_clustered);
    colormap(colors);
    colorbar;
    caxis([-1, 1]);
    
    title('Cerebrális Autoregulációs Metrikák Korrelációs Mátrixa', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Tengelyek címkézése (ha nem túl sok metrika)
    if length(names_clustered) <= 50
        set(gca, 'XTick', 1:length(names_clustered), 'XTickLabel', names_clustered, ...
            'XTickLabelRotation', 45);
        set(gca, 'YTick', 1:length(names_clustered), 'YTickLabel', names_clustered);
    else
        xlabel('Metrikák', 'FontSize', 12);
        ylabel('Metrikák', 'FontSize', 12);
    end
    
    % Grid
    grid on;
    set(gca, 'GridAlpha', 0.3);
    
    % Save
    savefig('correlation_heatmap.fig');
    print('correlation_heatmap.png', '-dpng', '-r300');
    fprintf('Heatmap mentve: correlation_heatmap.png\n');
end

function band_analysis = analyze_band_correlations(corr_matrix, metric_names)
    % Sáv-specifikus korrelációs analízis
    
    bands = {'Endothelial', 'Neurogenic', 'Myogenic', 'Respiratory', 'Cardiac', 'VLF', 'LF'};
    band_analysis = struct();
    
    for i = 1:length(bands)
        band_name = bands{i};
        
        % Sávhoz tartozó metrikák keresése
        band_mask = contains(metric_names, band_name);
        band_metrics = metric_names(band_mask);
        
        if sum(band_mask) >= 2
            % Sávon belüli korrelációk
            band_corr = corr_matrix(band_mask, band_mask);
            
            % Átlagos intra-band korreláció
            upper_tri = triu(band_corr, 1);
            valid_corr = upper_tri(upper_tri ~= 0);
            
            band_analysis.(band_name).metrics = band_metrics;
            band_analysis.(band_name).mean_correlation = mean(abs(valid_corr));
            band_analysis.(band_name).max_correlation = max(abs(valid_corr));
            band_analysis.(band_name).correlation_matrix = band_corr;
            
            fprintf('  %s sáv: %d metrika, átlag |r| = %.3f\n', ...
                band_name, length(band_metrics), band_analysis.(band_name).mean_correlation);
        end
    end
end

function autoreg_analysis = analyze_autoregulation_correlations(corr_matrix, metric_names)
    % Autoregulációs metrikák korrelációs analízise
    
    % Autoregulációs kulcsszavak
    autoreg_keywords = {'COx', 'autoregulation', 'Bilateral', 'CVRI', 'Mx', 'PRx', 'Sx', ...
                        'HRx', 'TOHRx', 'Autoregulation', 'HAI'};
    
    % Autoregulációs metrikák keresése
    autoreg_mask = false(size(metric_names));
    for i = 1:length(autoreg_keywords)
        autoreg_mask = autoreg_mask | contains(metric_names, autoreg_keywords{i}, 'IgnoreCase', true);
    end
    
    autoreg_metrics = metric_names(autoreg_mask);
    autoreg_corr = corr_matrix(autoreg_mask, autoreg_mask);
    
    % Elemzés
    upper_tri = triu(autoreg_corr, 1);
    valid_corr = upper_tri(upper_tri ~= 0);
    
    autoreg_analysis.metrics = autoreg_metrics;
    autoreg_analysis.correlation_matrix = autoreg_corr;
    autoreg_analysis.mean_correlation = mean(abs(valid_corr));
    autoreg_analysis.std_correlation = std(abs(valid_corr));
    autoreg_analysis.num_strong_correlations = sum(abs(valid_corr) > 0.7);
    
    fprintf('Autoregulációs metrikák: %d\n', length(autoreg_metrics));
    fprintf('  Átlagos |korreláció|: %.3f ± %.3f\n', ...
        autoreg_analysis.mean_correlation, autoreg_analysis.std_correlation);
    fprintf('  Erős korrelációk (|r| > 0.7): %d\n', autoreg_analysis.num_strong_correlations);
end

function generate_correlation_report(results)
    % Részletes korreláció jelentés generálása
    
    filename = 'correlation_analysis_report.txt';
    fid = fopen(filename, 'w');
    
    fprintf(fid, '📊 CEREBRÁLIS AUTOREGULÁCIÓS METRIKÁK KORRELÁCIÓS ELEMZÉS\n');
    fprintf(fid, '========================================================\n\n');
    fprintf(fid, 'Generálva: %s\n', datestr(now));
    fprintf(fid, 'Elemzett metrikák száma: %d\n', length(results.metric_names));
    fprintf(fid, 'Páciensek száma: %d\n\n', size(results.data_clean, 1));
    
    % Top korrelációk
    fprintf(fid, '🔝 TOP 20 LEGERŐSEBB KORRELÁCIÓ\n');
    fprintf(fid, '===============================\n');
    top_correlations = results.top_correlations;
    for i = 1:min(20, height(top_correlations))
        fprintf(fid, '%2d. %s ↔ %s: r = %+.4f\n', i, ...
            top_correlations.Metric1{i}, top_correlations.Metric2{i}, ...
            top_correlations.Correlation(i));
    end
    fprintf(fid, '\n');
    
    % Sáv-specifikus elemzés
    fprintf(fid, '🌊 SÁVI KORRELÁCIÓS ELEMZÉS\n');
    fprintf(fid, '===========================\n');
    bands = fieldnames(results.band_analysis);
    for i = 1:length(bands)
        band_name = bands{i};
        band_data = results.band_analysis.(band_name);
        fprintf(fid, '%s SÁV:\n', upper(band_name));
        fprintf(fid, '  Metrikák száma: %d\n', length(band_data.metrics));
        fprintf(fid, '  Átlagos |korreláció|: %.3f\n', band_data.mean_correlation);
        fprintf(fid, '  Maximális |korreláció|: %.3f\n', band_data.max_correlation);
        fprintf(fid, '\n');
    end
    
    % Autoregulációs metrikák elemzése
    fprintf(fid, '🧠 AUTOREGULÁCIÓS METRIKÁK KORRELÁCIÓJA\n');
    fprintf(fid, '======================================\n');
    autoreg = results.autoregulation_analysis;
    fprintf(fid, 'Autoregulációs metrikák száma: %d\n', length(autoreg.metrics));
    fprintf(fid, 'Átlagos |korreláció|: %.3f ± %.3f\n', ...
        autoreg.mean_correlation, autoreg.std_correlation);
    fprintf(fid, 'Erős korrelációk (|r| > 0.7): %d\n', autoreg.num_strong_correlations);
    fprintf(fid, 'Közepes korrelációk (0.4 < |r| < 0.7): %d\n', ...
        sum(abs(autoreg.correlation_matrix(:)) > 0.4 & abs(autoreg.correlation_matrix(:)) <= 0.7)/2);
    fprintf(fid, '\n');
    
    % Metrika kategória elemzés
    fprintf(fid, '📈 METRIKA KATEGÓRIÁK KORRELÁCIÓS MINTÁZATAI\n');
    fprintf(fid, '===========================================\n');
    
    % COx családú metrikák
    cox_keywords = {'COx', 'Mx', 'PRx', 'Sx'};
    cox_mask = false(size(results.metric_names));
    for k = 1:length(cox_keywords)
        cox_mask = cox_mask | contains(results.metric_names, cox_keywords{k});
    end
    
    if sum(cox_mask) > 1
        cox_corr_matrix = results.correlation_matrix(cox_mask, cox_mask);
        cox_upper_tri = triu(cox_corr_matrix, 1);
        cox_correlations = cox_upper_tri(cox_upper_tri ~= 0);
        fprintf(fid, 'COx CSALÁDÚ METRIKÁK:\n');
        fprintf(fid, '  Metrikák száma: %d\n', sum(cox_mask));
        fprintf(fid, '  Átlagos inter-korreláció: %.3f ± %.3f\n', ...
            mean(abs(cox_correlations)), std(abs(cox_correlations)));
        fprintf(fid, '  Erős belső korrelációk: %d\n', sum(abs(cox_correlations) > 0.7));
        fprintf(fid, '\n');
    end
    
    % Komplexitási metrikák
    complexity_keywords = {'DFA', 'Hurst', 'MSE', 'Entropy', 'Fractal', 'Complexity'};
    complexity_mask = false(size(results.metric_names));
    for k = 1:length(complexity_keywords)
        complexity_mask = complexity_mask | contains(results.metric_names, complexity_keywords{k});
    end
    
    if sum(complexity_mask) > 1
        comp_corr_matrix = results.correlation_matrix(complexity_mask, complexity_mask);
        comp_upper_tri = triu(comp_corr_matrix, 1);
        comp_correlations = comp_upper_tri(comp_upper_tri ~= 0);
        fprintf(fid, 'KOMPLEXITÁSI METRIKÁK:\n');
        fprintf(fid, '  Metrikák száma: %d\n', sum(complexity_mask));
        fprintf(fid, '  Átlagos inter-korreláció: %.3f ± %.3f\n', ...
            mean(abs(comp_correlations)), std(abs(comp_correlations)));
        fprintf(fid, '  Erős belső korrelációk: %d\n', sum(abs(comp_correlations) > 0.7));
        fprintf(fid, '\n');
    end
    
    % Ajánlások
    fprintf(fid, '💡 KLINIKAI AJÁNLÁSOK\n');
    fprintf(fid, '=====================\n');
    fprintf(fid, '1. REDUNDÁNS METRIKÁK: Erős korrelációjú (r > 0.9) metrikák közül\n');
    fprintf(fid, '   válasszon ki egy reprezentatívat a klinikai használatra.\n\n');
    
    fprintf(fid, '2. KOMPLEMENTER METRIKÁK: Gyenge korrelációjú (|r| < 0.3) metrikák\n');
    fprintf(fid, '   különböző aspektusokat mérnek - kombinált használat javasolt.\n\n');
    
    fprintf(fid, '3. VALIDÁCIÓ: Közepes korrelációjú (0.4 < |r| < 0.7) metrikák\n');
    fprintf(fid, '   kölcsönösen validálhatják egymást.\n\n');
    
    fprintf(fid, '4. SÁVI SPECIFICITÁS: Azonos sávú metrikák várhatóan korrelálnak,\n');
    fprintf(fid, '   különböző sávú metrikák gyengébb korrelációja normális.\n\n');
    
    fclose(fid);
    fprintf('Részletes korreláció jelentés: %s\n', filename);
end

function perform_advanced_correlation_analysis(results)
    % Haladó korrelációs elemzések
    
    fprintf('\n=== HALADÓ KORRELÁCIÓS ELEMZÉSEK ===\n');
    
    % 1. Principal Component Analysis
    fprintf('Principal Component Analysis végrehajtása...\n');
    try
        [coeff, score, latent, tsquared, explained] = pca(results.data_clean, 'Rows', 'pairwise');
        
        % PC eredmények mentése
        pc_results = struct();
        pc_results.coefficients = coeff;
        pc_results.scores = score;
        pc_results.explained_variance = explained;
        pc_results.cumsum_explained = cumsum(explained);
        
        % PC plot generálása
        figure('Position', [200, 200, 1000, 600]);
        
        subplot(1, 2, 1);
        bar(explained(1:20));
        title('Explained Variance by Principal Components (Top 20)');
        xlabel('Principal Component');
        ylabel('Explained Variance (%)');
        grid on;
        
        subplot(1, 2, 2);
        plot(cumsum(explained), 'LineWidth', 2);
        title('Cumulative Explained Variance');
        xlabel('Number of Components');
        ylabel('Cumulative Explained Variance (%)');
        grid on;
        yline(80, '--r', '80% Threshold');
        yline(90, '--g', '90% Threshold');
        
        savefig('pca_analysis.fig');
        print('pca_analysis.png', '-dpng', '-r300');
        
        % PC80 és PC90 számítása
        pc80_count = find(cumsum(explained) >= 80, 1);
        pc90_count = find(cumsum(explained) >= 90, 1);
        
        fprintf('  80%% variancia: %d komponens\n', pc80_count);
        fprintf('  90%% variancia: %d komponens\n', pc90_count);
        
        results.pca_analysis = pc_results;
        
    catch ME
        fprintf('  PCA hiba: %s\n', ME.message);
    end
    
    % 2. Hierarchikus klaszterezés
    fprintf('Hierarchikus klaszterezés végrehajtása...\n');
    try
        % Távolság mátrix (1 - |correlation|)
        dist_matrix = 1 - abs(results.correlation_matrix);
        
        % Linkage számítás
        linkage_matrix = linkage(squareform(dist_matrix), 'average');
        
        % Dendrogram
        figure('Position', [300, 300, 1200, 800]);
        
        if length(results.metric_names) <= 50
            % Ha kevés metrika, mutassuk a neveket
            dendrogram(linkage_matrix, 'Labels', results.metric_names, ...
                'Orientation', 'left');
        else
            % Ha sok metrika, csak a struktura
            dendrogram(linkage_matrix, 50);
            title('Hierarchical Clustering of Metrics (Top 50 branches)');
        end
        
        savefig('hierarchical_clustering.fig');
        print('hierarchical_clustering.png', '-dpng', '-r300');
        
        % Klaszterek képzése (pl. 10 klaszter)
        clusters = cluster(linkage_matrix, 'MaxClust', 10);
        results.clusters = clusters;
        
        fprintf('  10 klaszter létrehozva\n');
        
    catch ME
        fprintf('  Klaszterezés hiba: %s\n', ME.message);
    end
    
    % 3. Network analízis
    fprintf('Network analízis végrehajtása...\n');
    try
        % Strong correlation network (|r| > 0.5)
        strong_corr_threshold = 0.5;
        adj_matrix = abs(results.correlation_matrix) > strong_corr_threshold;
        
        % Remove self-connections
        adj_matrix = adj_matrix - diag(diag(adj_matrix));
        
        % Network properties
        network_stats = struct();
        network_stats.num_edges = sum(adj_matrix(:)) / 2;
        network_stats.density = network_stats.num_edges / (length(results.metric_names) * (length(results.metric_names) - 1) / 2);
        network_stats.degree_distribution = sum(adj_matrix, 2);
        network_stats.max_degree = max(network_stats.degree_distribution);
        network_stats.mean_degree = mean(network_stats.degree_distribution);
        
        fprintf('  Network edges (|r| > %.1f): %d\n', strong_corr_threshold, network_stats.num_edges);
        fprintf('  Network density: %.3f\n', network_stats.density);
        fprintf('  Average degree: %.1f\n', network_stats.mean_degree);
        fprintf('  Max degree: %d\n', network_stats.max_degree);
        
        % Hub metrikák (legnagyobb degree)
        [~, hub_indices] = sort(network_stats.degree_distribution, 'descend');
        fprintf('  Top 5 hub metrikák:\n');
        for i = 1:min(5, length(hub_indices))
            idx = hub_indices(i);
            fprintf('    %s (degree: %d)\n', results.metric_names{idx}, network_stats.degree_distribution(idx));
        end
        
        results.network_analysis = network_stats;
        
    catch ME
        fprintf('  Network analízis hiba: %s\n', ME.message);
    end
    
    % Eredmények mentése
    save('advanced_correlation_analysis.mat', 'results');
    fprintf('Haladó elemzés eredményei mentve: advanced_correlation_analysis.mat\n');
end

% Main függvény kiegészítése - add hozzá a main function végére
function correlation_analysis_results = analyze_metric_correlations()
    fprintf('=== CEREBRÁLIS AUTOREGULÁCIÓS METRIKÁK KORRELÁCIÓS ANALÍZISE ===\n\n');
    
    % CSV beolvasása
    try
        data_table = readtable('specific_indices_matrix.csv');
        fprintf('CSV sikeresen beolvasva: %d páciensek × %d változók\n', height(data_table), width(data_table));
    catch
        error('Hiba a CSV beolvasásakor. Ellenőrizd hogy a specific_indices_matrix.csv létezik!');
    end
    
    % PatientID oszlop eltávolítása
    if ismember('PatientID', data_table.Properties.VariableNames)
        data_numeric = data_table(:, 2:end); % PatientID kihagyása
        fprintf('PatientID oszlop eltávolítva. Elemzendő metrikák: %d\n', width(data_numeric));
    else
        data_numeric = data_table;
    end
    
    % Numerikus mátrix konvertálás
    data_matrix = table2array(data_numeric);
    metric_names = data_numeric.Properties.VariableNames;
    
    % NaN értékek ellenőrzése
    nan_counts = sum(isnan(data_matrix), 1);
    valid_metrics = nan_counts < height(data_table) * 0.8; % Max 80% NaN
    
    fprintf('Metrikák NaN eloszlása:\n');
    fprintf('  Teljes metrikák: %d\n', sum(valid_metrics));
    fprintf('  Túl sok NaN (>80%%): %d\n', sum(~valid_metrics));
    
    % Csak valid metrikák megtartása
    data_clean = data_matrix(:, valid_metrics);
    metrics_clean = metric_names(valid_metrics);
    
    fprintf('\nFinális analízis: %d metrika\n', length(metrics_clean));
    
    % Korrelációs mátrix számítása
    fprintf('\nKorrelációs mátrix számítása...\n');
    correlation_matrix = corr(data_clean, 'rows', 'pairwise');
    
    % Eredmények struct
    results = struct();
    results.data_table = data_table;
    results.correlation_matrix = correlation_matrix;
    results.metric_names = metrics_clean;
    results.data_clean = data_clean;
    
    % Top korrelációk keresése
    fprintf('Top korrelációk elemzése...\n');
    top_correlations = find_top_correlations(correlation_matrix, metrics_clean);
    results.top_correlations = top_correlations;
    
    % Heatmap generálása
    fprintf('Korrelációs heatmap generálása...\n');
    generate_correlation_heatmap(correlation_matrix, metrics_clean);
    
    % Sáv-specifikus korrelációs analízis
    fprintf('Sáv-specifikus korrelációs analízis...\n');
    band_analysis = analyze_band_correlations(correlation_matrix, metrics_clean);
    results.band_analysis = band_analysis;
    
    % Autoregulációs metrikák korrelációja
    fprintf('Autoregulációs metrikák korrelációja...\n');
    autoreg_analysis = analyze_autoregulation_correlations(correlation_matrix, metrics_clean);
    results.autoregulation_analysis = autoreg_analysis;
    
    % Részletes jelentés generálása
    fprintf('Részletes jelentés generálása...\n');
    generate_correlation_report(results);
    
    % Haladó elemzések
    fprintf('Haladó elemzések végrehajtása...\n');
    perform_advanced_correlation_analysis(results);
    
    % Eredmények mentése
    save('correlation_analysis_results.mat', 'results');
    fprintf('\n=== KORRELÁCIÓS ANALÍZIS BEFEJEZVE ===\n');
    fprintf('Eredmények mentve: correlation_analysis_results.mat\n');
    fprintf('Generált fájlok:\n');
    fprintf('  - correlation_heatmap.png\n');
    fprintf('  - top_correlations.csv\n');
    fprintf('  - correlation_analysis_report.txt\n');
    fprintf('  - pca_analysis.png\n');
    fprintf('  - hierarchical_clustering.png\n');
    fprintf('  - advanced_correlation_analysis.mat\n');
    
    correlation_analysis_results = results;
end
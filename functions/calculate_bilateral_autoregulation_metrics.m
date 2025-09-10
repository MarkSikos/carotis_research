function bilateral_results = calculate_bilateral_autoregulation_metrics(RSO2_ips, RSO2_contra, MAP_signal, fs)
    % Bilateral és laterality autoregulációs indexek
    % RSO2_ips = ipsilaterális (operált oldal), RSO2_contra = kontralaterális (másik oldal)
    
    bilateral_results = struct();
    fprintf('DEBUG: RSO2_ips std=%.3f, RSO2_contra std=%.3f, MAP std=%.3f\n', std(RSO2_ips), std(RSO2_contra), std(MAP_signal));
    fprintf('DEBUG: calculate_COx_single_hemisphere exists: %d\n', exist('calculate_COx_single_hemisphere', 'file'));
    
    try
        %% 1.1 BILATERAL COx (bCOx) - VALÓDI BILATERAL ADATOKKAL
        % Ipsilateral (operált oldal) COx
        COx_ipsi = calculate_COx_single_hemisphere(RSO2_ips, MAP_signal, fs);
        % Contralateral (másik oldal) COx
        COx_contra = calculate_COx_single_hemisphere(RSO2_contra, MAP_signal, fs);
        
        bilateral_results.COx_left = COx_ipsi;      % Ipsilateral mint "left"
        bilateral_results.COx_right = COx_contra;   % Contralateral mint "right"
        bilateral_results.COx_bilateral_mean = (COx_ipsi + COx_contra) / 2;
        bilateral_results.COx_bilateral_diff = abs(COx_ipsi - COx_contra);
        
        %% 1.2 HEMISPHERIC ASYMMETRY INDEX (HAI)
        % Aszimmetria index >10% klinikailag szignifikáns
        if abs(COx_ipsi) + abs(COx_contra) > 0
            bilateral_results.HAI_COx = abs(COx_ipsi - COx_contra) / ((abs(COx_ipsi) + abs(COx_contra))/2) * 100;
        else
            bilateral_results.HAI_COx = 0;
        end
        
        %% 1.3 IPSILATERAL vs CONTRALATERAL ANALÍZIS
        bilateral_results.COx_ipsilateral = COx_ipsi;
        bilateral_results.COx_contralateral = COx_contra;
        bilateral_results.laterality_index = COx_ipsi - COx_contra; % Pozitív = ipsi rosszabb
        
        %% 1.4 ASYMMETRIC BLOOD FLOW INDEX (ASYMrBF) - VALÓDI BILATERAL ADATOKKAL
        % Bilateral rSO2 flow különbség quantification
        if std(RSO2_ips) > 0 && std(RSO2_contra) > 0
            cross_corr_bilateral = corrcoef(RSO2_ips, RSO2_contra);
            bilateral_results.ASYMrBF = 1 - abs(cross_corr_bilateral(1,2)); % 0 = tökéletes szimmetria
        else
            bilateral_results.ASYMrBF = NaN;
        end
        
        %% 1.5 BILATERAL AUTOREGULATION EFFICIENCY
        % Kombinált bilateral autoregulációs hatékonyság
        if isfinite(COx_ipsi) && isfinite(COx_contra)
            bilateral_efficiency = (1 - abs(COx_ipsi)) * (1 - abs(COx_contra));
            bilateral_results.bilateral_autoregulation_efficiency = bilateral_efficiency;
        else
            bilateral_results.bilateral_autoregulation_efficiency = NaN;
        end
        
        fprintf('DEBUG: COx_ipsi=%.3f, COx_contra=%.3f, HAI=%.3f, Efficiency=%.3f\n', ...
                COx_ipsi, COx_contra, bilateral_results.HAI_COx, bilateral_results.bilateral_autoregulation_efficiency);
        
    catch ME
        fprintf('ERROR in bilateral metrics: %s\n', ME.message);
        % Error handling
        fields = {'COx_left', 'COx_right', 'COx_bilateral_mean', 'COx_bilateral_diff', ...
                  'HAI_COx', 'COx_ipsilateral', 'COx_contralateral', 'laterality_index', ...
                  'ASYMrBF', 'bilateral_autoregulation_efficiency'};
        for i = 1:length(fields)
            bilateral_results.(fields{i}) = NaN;
        end
    end
end
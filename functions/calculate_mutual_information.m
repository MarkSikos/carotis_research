function mi = calculate_mutual_information(signal1, signal2)
    % Mutual Information számítás
    try
        % Kvantálás (binning)
        n_bins = min(20, round(sqrt(length(signal1))));
        
        % Hisztogram edges
        edges1 = linspace(min(signal1), max(signal1), n_bins+1);
        edges2 = linspace(min(signal2), max(signal2), n_bins+1);
        
        % 2D hisztogram
        [N, ~, ~] = histcounts2(signal1, signal2, edges1, edges2);
        
        % Normalizálás
        P_xy = N / sum(N(:));
        P_x = sum(P_xy, 2);
        P_y = sum(P_xy, 1);
        
        % Mutual information számítás
        mi = 0;
        for i = 1:length(P_x)
            for j = 1:length(P_y)
                if P_xy(i,j) > 0 && P_x(i) > 0 && P_y(j) > 0
                    mi = mi + P_xy(i,j) * log2(P_xy(i,j) / (P_x(i) * P_y(j)));
                end
            end
        end
        
    catch
        mi = NaN;
    end
end


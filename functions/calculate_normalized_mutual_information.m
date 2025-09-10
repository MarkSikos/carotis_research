function nmi = calculate_normalized_mutual_information(signal1, signal2)
    % Normalized Mutual Information
    try
        mi = calculate_mutual_information(signal1, signal2);
        
        % Entropy calculation
        h1 = calculate_entropy(signal1);
        h2 = calculate_entropy(signal2);
        
        if h1 > 0 && h2 > 0
            nmi = 2 * mi / (h1 + h2);
        else
            nmi = NaN;
        end
        
    catch
        nmi = NaN;
    end
    
    function entropy = calculate_entropy(signal)
        n_bins = min(20, round(sqrt(length(signal))));
        [counts, ~] = histcounts(signal, n_bins);
        p = counts / sum(counts);
        p = p(p > 0); % Remove zero probabilities
        entropy = -sum(p .* log2(p));
    end
end

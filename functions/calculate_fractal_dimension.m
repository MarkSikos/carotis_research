function fractal_dim = calculate_fractal_dimension(signal)
    % Higuchi Fractal Dimension
    try
        signal = signal(:);
        N = length(signal);
        
        if N < 20
            fractal_dim = NaN;
            return;
        end
        
        % Higuchi algoritmus
        kmax = min(10, floor(N/4));
        k_values = 1:kmax;
        L_k = zeros(length(k_values), 1);
        
        for idx = 1:length(k_values)
            k = k_values(idx);
            L_m = zeros(k, 1);
            
            for m = 1:k
                indices = m:k:N;
                if length(indices) > 1
                    X_m = signal(indices);
                    L_m(m) = sum(abs(diff(X_m))) * (N - 1) / (length(indices) - 1) / k;
                end
            end
            
            L_k(idx) = mean(L_m(L_m > 0));
        end
        
        % Linear fit a log-log plot-ban
        valid_idx = L_k > 0 & isfinite(L_k);
        if sum(valid_idx) >= 3
            log_k = log10(k_values(valid_idx)');
            log_L = log10(L_k(valid_idx));
            
            p = polyfit(log_k, log_L, 1);
            fractal_dim = -p(1); % Negatív meredekség
        else
            fractal_dim = NaN;
        end
        
    catch
        fractal_dim = NaN;
    end
end

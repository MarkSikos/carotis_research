
function hurst = calculate_hurst_exponent(signal)
    % Hurst kitevo számítás R/S analízissel
    try
        signal = signal(:);
        N = length(signal);
        
        if N < 20
            hurst = NaN;
            return;
        end
        
        % R/S analízis
        max_tau = floor(N/4);
        tau_values = 2:max_tau;
        RS = zeros(length(tau_values), 1);
        
        for i = 1:length(tau_values)
            tau = tau_values(i);
            n_segments = floor(N / tau);
            rs_segments = zeros(n_segments, 1);
            
            for j = 1:n_segments
                start_idx = (j-1) * tau + 1;
                end_idx = j * tau;
                segment = signal(start_idx:end_idx);
                
                % Mean-centered cumulative sum
                mean_segment = mean(segment);
                Y = cumsum(segment - mean_segment);
                
                % Range
                R = max(Y) - min(Y);
                
                % Standard deviation  
                S = std(segment);
                
                if S > 0
                    rs_segments(j) = R / S;
                end
            end
            
            RS(i) = mean(rs_segments(rs_segments > 0));
        end
        
        % Linear fit a log-log plot-ban
        valid_idx = RS > 0 & isfinite(RS);
        if sum(valid_idx) >= 3
            log_tau = log10(tau_values(valid_idx)');
            log_RS = log10(RS(valid_idx));
            
            p = polyfit(log_tau, log_RS, 1);
            hurst = p(1); % Meredekség = Hurst kitevő
        else
            hurst = NaN;
        end
        
    catch
        hurst = NaN;
    end
end


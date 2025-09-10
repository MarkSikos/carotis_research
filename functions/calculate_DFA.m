function results = calculate_DFA(SE_signal, rSO2_signal)
    % Detrended Fluctuation Analysis
    results = struct();
    try
        results.SE_alpha = DFA_single(SE_signal);
        results.rSO2_alpha = DFA_single(rSO2_signal);
        % Cross-DFA (simplified)
        results.cross_DFA = (results.SE_alpha + results.rSO2_alpha) / 2;
    catch
        results.SE_alpha = NaN;
        results.rSO2_alpha = NaN;
        results.cross_DFA = NaN;
    end

    function alpha = DFA_single(signal)
        signal = signal(:);
        N = length(signal);
        if N < 20
            alpha = NaN;
            return;
        end
        % Integration
        y = cumsum(signal - mean(signal));
        % Box sizes
        min_box = 4;
        max_box = floor(N/4);
        box_sizes = round(logspace(log10(min_box), log10(max_box), 15));
        box_sizes = unique(box_sizes);
        F = zeros(length(box_sizes), 1);
        for i = 1:length(box_sizes)
            n = box_sizes(i);
            n_boxes = floor(N / n);
            mse = 0;
            for j = 1:n_boxes
                start_idx = (j-1) * n + 1;
                end_idx = j * n;
                segment = y(start_idx:end_idx);
                x = (1:n)';
                % Linear detrending
                p = polyfit(x, segment, 1);
                trend = polyval(p, x);
                mse = mse + mean((segment - trend).^2);
            end
            F(i) = sqrt(mse / n_boxes);
        end
        % Linear fit a log-log plot-ban
        valid_idx = F > 0 & isfinite(F);
        if sum(valid_idx) >= 3
            log_n = log10(box_sizes(valid_idx)');
            log_F = log10(F(valid_idx));
            p = polyfit(log_n, log_F, 1);
            alpha = p(1);
        else
            alpha = NaN;
        end
    end
end
function results = calculate_multiscale_entropy(SE_signal, rSO2_signal)
    % Multiscale Entropy számítás
    results = struct();
    try
        max_scale = min(20, floor(length(SE_signal)/10));
        scales = 1:max_scale;
        SE_MSE = zeros(length(scales), 1);
        rSO2_MSE = zeros(length(scales), 1);
        for i = 1:length(scales)
            scale = scales(i);
            % Coarse-graining
            SE_coarse = coarse_grain(SE_signal, scale);
            rSO2_coarse = coarse_grain(rSO2_signal, scale);
            % Sample entropy az adott skálán
            SE_MSE(i) = calculate_sample_entropy(SE_coarse);
            rSO2_MSE(i) = calculate_sample_entropy(rSO2_coarse);
        end
        results.SE_MSE = mean(SE_MSE(isfinite(SE_MSE)));
        results.rSO2_MSE = mean(rSO2_MSE(isfinite(rSO2_MSE)));
        % Complexity Index (terület a MSE görbe alatt)
        valid_SE = SE_MSE(isfinite(SE_MSE));
        valid_rSO2 = rSO2_MSE(isfinite(rSO2_MSE));
        if length(valid_SE) > 1 && length(valid_rSO2) > 1
            results.complexity_index = trapz(valid_SE) + trapz(valid_rSO2);
        else
            results.complexity_index = NaN;
        end
    catch
        results.SE_MSE = NaN;
        results.rSO2_MSE = NaN;
        results.complexity_index = NaN;
    end

    function coarse = coarse_grain(signal, scale)
        N = length(signal);
        coarse_length = floor(N / scale);
        coarse = zeros(coarse_length, 1);
        for i = 1:coarse_length
            start_idx = (i-1) * scale + 1;
            end_idx = i * scale;
            coarse(i) = mean(signal(start_idx:end_idx));
        end
    end

    function sampen = calculate_sample_entropy(signal)
        % Sample Entropy számítás
        try
            signal = signal(:);
            N = length(signal);
            if N < 10
                sampen = NaN;
                return;
            end
            
            m = 2; % Pattern length
            r = 0.2 * std(signal); % Tolerance
            
            if r == 0
                sampen = NaN;
                return;
            end
            
            % Count matches for m and m+1
            matches_m = 0;
            matches_m1 = 0;
            
            for i = 1:(N-m)
                template_m = signal(i:i+m-1);
                template_m1 = signal(i:i+m);
                
                for j = 1:(N-m)
                    if i ~= j
                        test_m = signal(j:j+m-1);
                        test_m1 = signal(j:j+m);
                        
                        % Check match for m
                        if max(abs(template_m - test_m)) <= r
                            matches_m = matches_m + 1;
                            
                            % Check match for m+1
                            if j <= N-m-1 && max(abs(template_m1 - test_m1)) <= r
                                matches_m1 = matches_m1 + 1;
                            end
                        end
                    end
                end
            end
            
            if matches_m > 0 && matches_m1 > 0
                sampen = -log(matches_m1 / matches_m);
            else
                sampen = NaN;
            end
            
        catch
            sampen = NaN;
        end
    end
end
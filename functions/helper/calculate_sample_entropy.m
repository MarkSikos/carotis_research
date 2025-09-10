function sampen = calculate_sample_entropy(signal)
    % Sample Entropy számítás
    try
        signal = signal(:);
        N = length(signal);
        
        if N < 10
            sampen = NaN;
            return;
        end
        
        % Paraméterek
        m = 2; % Pattern length
        r = 0.2 * std(signal); % Tolerance
        
        % Template matching
        function C = maxdist(xi, xj, N)
            C = max(abs(xi(1:N) - xj(1:N)));
        end
        
        phi = zeros(1, 2);
        
        for j = 1:2
            m_temp = m + j - 1;
            patterns = zeros(N - m_temp, m_temp + 1);
            
            for i = 1:(N - m_temp)
                patterns(i, :) = signal(i:i+m_temp);
            end
            
            C = 0;
            for i = 1:(N - m_temp)
                template = patterns(i, 1:m_temp);
                
                for k = 1:(N - m_temp)
                    if i ~= k
                        test_pattern = patterns(k, 1:m_temp);
                        if maxdist(template, test_pattern, m_temp) <= r
                            C = C + 1;
                        end
                    end
                end
            end
            
            phi(j) = C / ((N - m_temp) * (N - m_temp - 1));
        end
        
        if phi(1) > 0
            sampen = -log(phi(2) / phi(1));
        else
            sampen = NaN;
        end
        
    catch
        sampen = NaN;
    end
end

function results = calculate_surrogate_testing(SE_clean, rSO2_clean)
    % Surrogate data testing a kapcsolat szignifikanciájára
    results = struct();
    try
        % Eredeti korreláció
        original_corr = corr(SE_clean, rSO2_clean);
        % Surrogate adatok generálása
        n_surrogates = 100;
        surrogate_corrs = zeros(n_surrogates, 1);
        for i = 1:n_surrogates
            % Phase randomization (Fourier surrogate)
            rSO2_surrogate = generate_fourier_surrogate(rSO2_clean);
            surrogate_corrs(i) = corr(SE_clean, rSO2_surrogate);
        end
        % P-value számítás
        if abs(original_corr) > 0
            p_value = sum(abs(surrogate_corrs) >= abs(original_corr)) / n_surrogates;
            results.p_value = p_value;
            results.significance = p_value < 0.05;
        else
            results.p_value = 1;
            results.significance = false;
        end
    catch
        results.p_value = NaN;
        results.significance = false;
    end

    function surrogate = generate_fourier_surrogate(signal)
        % Phase randomization
        try
            fft_signal = fft(signal);
            magnitude = abs(fft_signal);
            random_phases = angle(fft_signal);
            random_phases(2:end-1) = 2*pi*rand(length(random_phases)-2, 1) - pi;
            % Hermitian symmetry preservation
            n = length(signal);
            if mod(n, 2) == 0
                random_phases(n/2+2:end) = -random_phases(n/2:-1:2);
            else
                random_phases((n+1)/2+1:end) = -random_phases((n+1)/2:-1:2);
            end
            surrogate_fft = magnitude .* exp(1i * random_phases);
            surrogate = real(ifft(surrogate_fft));
        catch
            surrogate = signal; % Fallback ha hiba van
        end
    end
end
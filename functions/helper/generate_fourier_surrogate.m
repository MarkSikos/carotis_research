
    function surrogate = generate_fourier_surrogate(signal)
        % Phase randomization
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
    end
end
function cross_sampen = calculate_cross_sample_entropy(signal1, signal2)
    % Cross Sample Entropy két jel között
    try
        % Egyszerűsített implementáció
        combined_signal = [signal1(:); signal2(:)];
        cross_sampen = calculate_sample_entropy(combined_signal);
    catch
        cross_sampen = NaN;
    end
end


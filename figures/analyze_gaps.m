

function gaps = analyze_gaps(signal)
% Elemzi a NaN gap-ek hosszát egy signalban
gaps = [];
i = 1;
while i <= length(signal)
    if isnan(signal(i))
        gap_start = i;
        while i <= length(signal) && isnan(signal(i))
            i = i + 1;
        end
        gap_length = i - gap_start;
        gaps = [gaps, gap_length];
    else
        i = i + 1;
    end
end
end
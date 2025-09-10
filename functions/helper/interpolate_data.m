function interpolated_data = interpolate_data(patient_data, n_points, target_duration)
    % Sáv-specifikus interpoláció SE, rSO2, MAP, HR, other_side_oxig és resistance jelekhez
    original_time = patient_data.Time;
    original_SE = patient_data.SE;
    original_rSO2 = patient_data.oper_side_oxig;
    original_MAP = patient_data.MAP;
    original_HR = patient_data.HR; % HR hozzáadása
    original_rSO2_contra = patient_data.other_side_oxig; % Contralateral rSO2 hozzáadása
    original_resistance = patient_data.resistance; % ÚJ: Resistance hozzáadása
    
    % Időintervallum kiválasztása
    total_duration = max(original_time) - min(original_time);
    if total_duration >= target_duration
        % Van elég adat, középső részt vesszük
        start_time = min(original_time) + (total_duration - target_duration) / 2;
        end_time = start_time + target_duration;
        time_mask = original_time >= start_time & original_time <= end_time;
        selected_time = original_time(time_mask);
        selected_SE = original_SE(time_mask);
        selected_rSO2 = original_rSO2(time_mask);
        selected_MAP = original_MAP(time_mask);
        selected_HR = original_HR(time_mask); % HR kiválasztása
        selected_rSO2_contra = original_rSO2_contra(time_mask); % Contralateral kiválasztása
        selected_resistance = original_resistance(time_mask); % ÚJ: Resistance kiválasztása
    else
        % Kevés adat, az egészet használjuk
        selected_time = original_time;
        selected_SE = original_SE;
        selected_rSO2 = original_rSO2;
        selected_MAP = original_MAP;
        selected_HR = original_HR; % HR használata
        selected_rSO2_contra = original_rSO2_contra; % Contralateral használata
        selected_resistance = original_resistance; % ÚJ: Resistance használata
    end
    
    % Interpoláció mind a hat jelre
    new_time = linspace(min(selected_time), max(selected_time), n_points);
    interpolated_data = struct();
    interpolated_data.SE = interp1(selected_time, selected_SE, new_time, 'pchip', 'extrap');
    interpolated_data.rSO2 = interp1(selected_time, selected_rSO2, new_time, 'pchip', 'extrap');
    interpolated_data.MAP = interp1(selected_time, selected_MAP, new_time, 'pchip', 'extrap');
    interpolated_data.HR = interp1(selected_time, selected_HR, new_time, 'pchip', 'extrap'); % HR interpoláció
    interpolated_data.rSO2_contra = interp1(selected_time, selected_rSO2_contra, new_time, 'pchip', 'extrap'); % Contralateral interpoláció
    interpolated_data.resistance = interp1(selected_time, selected_resistance, new_time, 'pchip', 'extrap'); % ÚJ: Resistance interpoláció
    interpolated_data.time = new_time;
end
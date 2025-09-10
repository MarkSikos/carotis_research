function results = calculate_granger_causality(SE_signal, rSO2_signal)
    % Granger Causality analízis (egyszerűsített)
    results = struct();
    try
        % Adatok előkészítése
        data = [SE_signal(:), rSO2_signal(:)];
        data = data(~any(isnan(data), 2), :); % NaN-ok eltávolítása
        if size(data, 1) < 20
            results.SE_to_rSO2 = NaN;
            results.rSO2_to_SE = NaN;
            results.bidirectional = NaN;
            return;
        end
        % VAR model order selection (egyszerű AIC alapú)
        max_order = min(10, floor(size(data, 1)/10));
        best_order = 1;
        best_aic = inf;
        for order = 1:max_order
            try
                % VAR(order) modell
                [~, aic] = estimate_var_model(data, order);
                if aic < best_aic
                    best_aic = aic;
                    best_order = order;
                end
            catch
                continue;
            end
        end
        % Granger teszt a legjobb order-rel
        order = best_order;
        % SE -> rSO2 irány
        gc_SE_to_rSO2 = granger_test_direction(data(:,1), data(:,2), order);
        % rSO2 -> SE irány
        gc_rSO2_to_SE = granger_test_direction(data(:,2), data(:,1), order);
        results.SE_to_rSO2 = gc_SE_to_rSO2;
        results.rSO2_to_SE = gc_rSO2_to_SE;
        results.bidirectional = (gc_SE_to_rSO2 + gc_rSO2_to_SE) / 2;
    catch
        results.SE_to_rSO2 = NaN;
        results.rSO2_to_SE = NaN;
        results.bidirectional = NaN;
    end

    function [model, aic] = estimate_var_model(data, order)
        % Egyszerűsített VAR model becslés
        try
            n = size(data, 1);
            k = size(data, 2);
            
            if n <= order * k + k
                aic = inf;
                model = [];
                return;
            end
            
            % Lagged variables építése
            Y = data(order+1:end, :);
            X = [];
            for lag = 1:order
                X = [X, data(order+1-lag:end-lag, :)];
            end
            X = [ones(size(Y, 1), 1), X]; % Konstans tag
            
            % OLS becslés
            beta = X \ Y;
            residuals = Y - X * beta;
            
            % AIC számítás
            rss = sum(sum(residuals.^2));
            n_params = size(X, 2) * k;
            aic = n * log(rss/n) + 2 * n_params;
            
            model.beta = beta;
            model.residuals = residuals;
        catch
            aic = inf;
            model = [];
        end
    end

    function gc_value = granger_test_direction(x, y, order)
        % Granger causality teszt: x -> y irány
        try
            data = [x(:), y(:)];
            data = data(~any(isnan(data), 2), :);
            
            if size(data, 1) <= order + 10
                gc_value = NaN;
                return;
            end
            
            x_clean = data(:, 1);
            y_clean = data(:, 2);
            n = length(y_clean);
            
            % Korlátozott modell: y(t) = c + sum(a_i * y(t-i))
            Y_restricted = y_clean(order+1:end);
            X_restricted = ones(length(Y_restricted), 1);
            for lag = 1:order
                X_restricted = [X_restricted, y_clean(order+1-lag:end-lag)];
            end
            
            beta_restricted = X_restricted \ Y_restricted;
            residuals_restricted = Y_restricted - X_restricted * beta_restricted;
            rss_restricted = sum(residuals_restricted.^2);
            
            % Teljes modell: y(t) = c + sum(a_i * y(t-i)) + sum(b_i * x(t-i))
            Y_full = y_clean(order+1:end);
            X_full = X_restricted;
            for lag = 1:order
                X_full = [X_full, x_clean(order+1-lag:end-lag)];
            end
            
            beta_full = X_full \ Y_full;
            residuals_full = Y_full - X_full * beta_full;
            rss_full = sum(residuals_full.^2);
            
            % F-teszt statisztika proxy (egyszerűsített)
            if rss_full > 0 && rss_restricted > rss_full
                gc_value = log(rss_restricted / rss_full);
            else
                gc_value = 0;
            end
            
        catch
            gc_value = NaN;
        end
    end
end
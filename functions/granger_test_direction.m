 function gc_value = granger_test_direction(cause, effect, order)
        try
            n = length(cause);
            
            % Restricted model: effect ~ effect_lags
            X_restricted = [];
            Y_restricted = [];
            
            for t = order+1:n
                y_t = effect(t);
                x_t = [];
                for lag = 1:order
                    x_t = [x_t; effect(t-lag)];
                end
                X_restricted = [X_restricted, x_t];
                Y_restricted = [Y_restricted, y_t];
            end
            
            % Unrestricted model: effect ~ effect_lags + cause_lags
            X_unrestricted = [];
            Y_unrestricted = [];
            
            for t = order+1:n
                y_t = effect(t);
                x_t = [];
                for lag = 1:order
                    x_t = [x_t; effect(t-lag); cause(t-lag)];
                end
                X_unrestricted = [X_unrestricted, x_t];
                Y_unrestricted = [Y_unrestricted, y_t];
            end
            
            % Residual sum of squares
            if size(X_restricted, 1) > 0 && size(X_unrestricted, 1) > 0
                beta_restricted = Y_restricted / X_restricted;
                beta_unrestricted = Y_unrestricted / X_unrestricted;
                
                RSS_restricted = sum((Y_restricted - beta_restricted * X_restricted).^2);
                RSS_unrestricted = sum((Y_unrestricted - beta_unrestricted * X_unrestricted).^2);
                
                % F-statistic approximáció
                if RSS_restricted > RSS_unrestricted && RSS_unrestricted > 0
                    gc_value = log(RSS_restricted / RSS_unrestricted);
                else
                    gc_value = 0;
                end
            else
                gc_value = NaN;
            end
            
        catch
            gc_value = NaN;
        end
    end
end


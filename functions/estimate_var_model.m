   function [coeffs, aic] = estimate_var_model(data, order)
        n = size(data, 1);
        p = size(data, 2);
        
        % Design matrix
        X = [];
        Y = [];
        
        for t = order+1:n
            y_t = data(t, :)';
            x_t = [];
            for lag = 1:order
                x_t = [x_t; data(t-lag, :)'];
            end
            X = [X, x_t];
            Y = [Y, y_t];
        end
        
        % OLS becslés
        coeffs = Y / X;
        
        % AIC számítás
        residuals = Y - coeffs * X;
        sigma = residuals * residuals' / size(residuals, 2);
        aic = log(det(sigma)) + 2 * order * p^2 / size(residuals, 2);
    end
    

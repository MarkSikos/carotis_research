
%% 6. MULTI-PARAMETER COMPOSITE AUTOREGULATION SCORE
%% ========================================================================

function composite_score = calculate_composite_autoregulation_score(COx, bilateral_results, resistance_results, HRx_results)
    % Composite autoregulation score combining multiple metrics
    
    try
        % Initialize score components
        scores = [];
        weights = [];
        
        % COx component (weight: 0.3)
        if isfinite(COx)
            cox_score = 1 - abs(COx); % Higher score = better autoregulation
            scores = [scores; cox_score];
            weights = [weights; 0.3];
        end
        
        % Bilateral component (weight: 0.25)
        if isfield(bilateral_results, 'bilateral_autoregulation_efficiency') && ...
           isfinite(bilateral_results.bilateral_autoregulation_efficiency)
            scores = [scores; bilateral_results.bilateral_autoregulation_efficiency];
            weights = [weights; 0.25];
        end
        
        % Resistance component (weight: 0.2)
        if isfield(resistance_results, 'CVRI') && isfinite(resistance_results.CVRI)
            resistance_score = 1 - abs(resistance_results.CVRI);
            scores = [scores; resistance_score];
            weights = [weights; 0.2];
        end
        
        % Heart rate component (weight: 0.15)
        if isfield(HRx_results, 'TOHRx') && isfinite(HRx_results.TOHRx)
            hr_score = 1 - abs(HRx_results.TOHRx);
            scores = [scores; hr_score];
            weights = [weights; 0.15];
        end
        
        % Asymmetry penalty (weight: 0.1)
        if isfield(bilateral_results, 'HAI_COx') && isfinite(bilateral_results.HAI_COx)
            asymmetry_penalty = bilateral_results.HAI_COx / 100; % Convert percentage to ratio
            asymmetry_score = 1 - min(asymmetry_penalty, 1); % Cap at 1
            scores = [scores; asymmetry_score];
            weights = [weights; 0.1];
        end
        
        % Calculate weighted composite score
        if ~isempty(scores)
            % Normalize weights to sum to 1
            weights = weights / sum(weights);
            composite_score = sum(scores .* weights);
            
            % Ensure score is between 0 and 1
            composite_score = max(0, min(1, composite_score));
        else
            composite_score = NaN;
        end
        
    catch
        composite_score = NaN;
    end
end
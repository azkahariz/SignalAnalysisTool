function errorThresh = findConvergenceThreshold(data, startVal, increment)
% FINDCONVERGENCETHRESHOLD Determines the error threshold for convergence.
%
% This function identifies the error threshold at which all Monte Carlo
% simulations within a SignalAnalysisTool object achieve convergence.
%
% Inputs:
%   - data: An object of type SignalAnalysisTool containing Monte Carlo 
%           simulation data and convergence criteria.
%   - startVal: The initial value of the error threshold to begin the search.
%   - increment: The step size for incrementing the error threshold during 
%                the search process.
%
% Outputs:
%   - errorThresh: The identified error threshold at which all simulations 
%                  achieve convergence.
%
% Example:
%   errorThresh = findConvergenceThreshold(data, 0.1, 0.01);
%
%   This example starts the search with an initial threshold of 0.1 and 
%   increments by 0.01 until all simulations converge.

    % Extract the total number of Monte Carlo iterations from the data object
    mcIter = data.MonteCarloIterations;
    
    % Initialize variables
    isConverged = false;  % Flag to indicate convergence of all simulations
    errorThresh = startVal;  % Start threshold value

    % Iteratively search for the convergence threshold
    while ~isConverged
        % Check how many simulations have converged at the current threshold
        numConverged = checkConvergence(data, errorThresh);
        
        % If all simulations are converged, terminate the search
        if numConverged == mcIter
            isConverged = true;
        else
            % Increment the threshold and continue
            errorThresh = errorThresh + increment;
        end
    end

    % Nested function to check the convergence of simulations
    function numConverged = checkConvergence(obj, thresh)
        % CHECKCONVERGENCE Checks the number of simulations converged.
        %
        % Inputs:
        %   - obj: The SignalAnalysisTool object containing simulation data.
        %   - thresh: The current error threshold to evaluate convergence.
        %
        % Outputs:
        %   - numConverged: The number of simulations that meet the 
        %                   convergence criteria.
        
        % Retrieve convergence status for x and y axes
        [xConv, yConv] = obj.checkErrorConvergence(thresh);
        
        % Count the number of simulations where both x and y are converged
        numConverged = sum(xConv & yConv);
    end

end
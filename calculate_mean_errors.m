function mean_errors = calculate_mean_errors(filter_data)
    % Function to calculate mean errors for a given filter data structure
    % INPUT:
    %   filter_data - Struct containing the filter's simulation data
    % OUTPUT:
    %   mean_errors - Row vector containing the mean errors for each state

    % Extract the true state from the simulation data
    x_true = filter_data.SimulationData.x.values{1};
    
    % Calculate the absolute error for each state
    errors = abs(filter_data.StatOfMonteCarlo.xMean - x_true);
    
    % Calculate the mean error for each state
    mean_errors = mean(errors, 1); % Mean along the rows
end

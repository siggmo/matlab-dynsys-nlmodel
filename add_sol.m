function sols = add_sol(sols, sol, T, u, h, title, color)
% ADD_SOL This function takes in a solution by an arbitrary matlab ode-solver and 
% adds it to an array suitable for plotting with plot_sols_2.m. Therefore one 
% needs to provide the existing array (sols), output map (h(x, u)), input
% function u(t), title of the dataset and the desired color for its plot,
% also the desired times T where the solution should be evaluated.
    
    % extract relevant data from ode-solvers result
    new_sol.x = deval(sol, T);  % states
    
    % evaluate at single points only, in case u or h don't support stacked
    % arguments. Determine dimensions and orientation first.
    h_temp = h(deval(sol, T(1)), u(T(1)));
    h_is_column = iscolumn(h_temp);

    % create zero arrays of suitable dimensions
    new_sol.y = zeros([length(h_temp) length(T)]);

    % evaluate h at each time step and store in previously created array
    for i=1:length(T)
        % make sure it's a column vector
        if h_is_column
            new_sol.y(:, i) = h(new_sol.x(:, i), u(T(i)));
        else
            new_sol.y(:, i) = h(new_sol.x(:, i), u(T(i)))';
        end
    end
    
    new_sol.T = T;
    new_sol.color = color;
    new_sol.desc = title;
    sols = [sols new_sol];  % add to sols
end
function plot_sols_2(plot_title, desc, select, limits, sols)
% Plot selected components of the solutions contained in
% sols.
% Descriptions and whether to show or not are expected to be
% found in desc and select arguments. All solutions contained in sols need
% to have the same number of components.

% sols is expected to be an array of objects of the following structure:
% sol.t: time vector
% sol.y: output data (components per row, samples per column)
% sol.color: color string for this dataset
% sol.desc: description of this dataset

% desc should contain description strings for each component.

% select is expected to be a 0/1 vector specifying for each component of 
% the solution whether it should be displayed or not.

    % analyze components
    n_components = length(desc);
    n_selected = sum(select);
    k = size(sols(1).y, 1);  % number of outputs in the solution
    
    % check dimensions
    if k ~= n_components
        disp('number of solution components and descriptions doesnt match!');
    end

    % plot
    figure('Name',plot_title,'NumberTitle','off');
    tiledlayout(n_selected+1, 1);
    %tiledlayout('flow');
    current_plot = 1;
    % loop through all components
    if n_components > 0
        if k > 0
            for i = 1:k
                if select(i)
                    nexttile(current_plot);
                    % loop through all solutions
                    % corresponding outputs of each solution are plotted in the
                    % same plot for comparison.
                    for j = 1:length(sols)
                        y_data = sols(j).y(i,:);
                        plot(sols(j).T, y_data, 'Color', sols(j).color, 'DisplayName', sols(j).desc);
                        hold on;
                    end
    
                    % only set limits if valid
                    if diff(limits(i, :)) > 0
                        ylim(limits(i, :));
                    end
                    
                    title(desc(i));hold on;grid on;
                    current_plot = current_plot + 1;
                end
            end
        end
    end

    % generate legend
    lgd = legend;
    lgd.Layout.Tile = current_plot;

end
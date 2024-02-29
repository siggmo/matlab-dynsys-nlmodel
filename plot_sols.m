function plot_sols(plot_title, output_desc, sols)
% simplified call of plot_sols_2. Shows all solution and input components
% and applies limits automatically.

m = length(output_desc);
show_outputs = ones([m 1]);
output_limits = zeros([m 2]);

plot_sols_2(plot_title, output_desc, show_outputs, output_limits, sols);
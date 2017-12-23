%% showing that (y-x) vs. (x) can have negative correlation
% if add random noise to y and x.


signal = 0.6;
noise = 0.2;

z = linspace(0, signal, 100);

x = noise*normrnd(0, noise, 1, 100) + z;
y = noise*normrnd(0, noise, 1, 100) + z;

lt_figure; hold on;
% ------------
lt_subplot(2,2,1); hold on;
title('ground truth values')
xlabel('x'); ylabel('y');

lt_plot_45degScatter(z,z, 'k')

% -----------
lt_subplot(2,2,2); hold on;
title('distribution of noise added to ground truth');
lt_plot_histogram(normrnd(0, noise, 1, 100));

% ------------
lt_subplot(2,2,3); hold on;
xlabel('x'); ylabel('y');

lt_regress(y,x,1);
lt_plot_45degScatter(x, y, 'k')
title('after noise added to x,y, for each datapoint');

% ------------
lt_subplot(2,2,4); hold on;
title('same (noised) data, but plotting y-x instead of y')
xlabel('x'); ylabel('y-x');

lt_regress(y-x, x, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
function [] = db_prepare_bootstrp_hist(fig_title,fig_ylabel,fig_xlabel,CI,rel_freq,median)
%db_preare_bootstrp_hist Prepares a figure for boostrap histogram data
%   Adds a title, x and y label, and lines (5th, 50th, and 95th percentile)
%   for histogram figure for boostrap data

title(fig_title,'interpreter','none')
ylabel(fig_ylabel)
xlabel(fig_xlabel)

% xlim([-.1 1.1]);
% ylim([-.1 1.1]);

% Line for 5th percentile
line([CI(1) CI(1)], [0 max(rel_freq)],...
    'LineWidth',3,'Color', [1 0 0])

% Line for 95th percentile
line([CI(2) CI(2)], [0 max(rel_freq)],...
    'LineWidth',3,'Color', [1 0 0])

% Line for 50th percentile (median)
line([median median], [0 max(rel_freq)],...
    'LineWidth',3,'Color', [0 1 0])



end


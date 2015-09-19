%% LT 8/20/15 - plots empirical cdf
function lt_plot_cdf(Y, color);


if ~exist('color','var');
    color='k';
end

%%
[F, X]=ecdf(Y);
plot(X, F, '-', 'Color',color,'LineWidth',2);


%% LT 8/20/15 - plots empirical cdf
function h=lt_plot_cdf(Y, color, plotCount);
if ~exist('plotCount', 'var');
    plotCount=0;
end

if ~exist('color','var');
    color='k';
end

%%
[F, X]=ecdf(Y);

if plotCount==1
 h=plot(X, F.*numel(Y), '-', 'Color',color,'LineWidth',2);
   
else
h=plot(X, F, '-', 'Color',color,'LineWidth',2);
end

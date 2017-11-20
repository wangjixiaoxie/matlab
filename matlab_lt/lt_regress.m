%% uses matlab code regress() but automatically adds 1 to X and formats X and Y. Also plots automatically
function [b,bint,r,rint,stats,SummaryStats, hplot]=lt_regress(Y,X_no_ones,plotON, ...
    plotB_error, plotline, plot_summary, plotcol, plotOnlyLineSummary, alpha)

% Y and X are vectors (can be columns or rows, currently only 1D) (do not
% add 1s to X)
% plotON =1 plots data with best fit line.
% plotB_error=1 plot 95% CI of line.
% plotline=1

% outputs: same as for regress()

%% DEFAULTS

if ~exist('plot_summary','var');
    plot_summary=1;
end

if ~exist('plotB_error','var');
    plotB_error=0;
end

if ~exist('plotline','var');
    plotline=1;
end

if ~exist('plotcol','var');
    plotcol='b';
end

if ~exist('plotOnlyLineSummary', 'var')
    plotOnlyLineSummary=0;
end

if ~exist('alpha', 'var');
    alpha = 0.05;
end

%% put inputs into columns
% assume they are 1d vectors

if size(Y,2)>size(Y,1);
    Y=Y';
end

if size(X_no_ones,2)>size(X_no_ones,1);
    X_no_ones=X_no_ones';
end

%% put ones in front of X
X=[ones(length(X_no_ones),1), X_no_ones];

%% perform regression

[b,bint,r,rint,stats]=regress(Y,X, alpha);

%% put some useful stuff into a summary structure

SummaryStats.slope=b(2);
SummaryStats.slopeCI = bint(2,:);
SummaryStats.intercept=b(1);
SummaryStats.interceptCI = bint(1,:);
SummaryStats.R2=stats(1);
SummaryStats.p=stats(3);


%% plot if desired
if plotON==1;
    %     lt_figure; hold on;
    %     title('raw data + linear regression best fit (95% CI)')
    % plot data
    hplot=plot(X(:,2),Y, '.', 'Color',plotcol);
    
    % plot regression lines (with error)
    if plotline==1;
        plot(xlim,b(1) + b(2).*xlim,'-','Color',plotcol,'LineWidth',2);
    end
    
    if plotB_error==1;
        plot(xlim,bint(1,1) + bint(2,1).*xlim,'-r');
        plot(xlim,bint(1,2) + bint(2,2).*xlim,'-r');
    end
    
    plotcolsummary = 'r';
    if plotcol =='r';
        plotcolsummary = 'k';
    end
    
    % plot summary stats
    if plot_summary==1;
        Xlim=xlim;
        Ylim=ylim;
        
        text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-(Ylim(2)-Ylim(1))/10,['p=' num2str(SummaryStats.p)],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
        text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-2*(Ylim(2)-Ylim(1))/10,['R2=' num2str(SummaryStats.R2)],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
        text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-3*(Ylim(2)-Ylim(1))/10,['slope=' num2str(SummaryStats.slope) '(' num2str(SummaryStats.slopeCI) ')'],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
        text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-4*(Ylim(2)-Ylim(1))/10,['int=' num2str(SummaryStats.intercept) '(' num2str(SummaryStats.interceptCI) ')'],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
        
    end
end

    plotcolsummary = 'r';
    if plotcol =='r';
        plotcolsummary = 'k';
    end

    
if plotOnlyLineSummary==1
    plot(xlim,b(1) + b(2).*xlim,'-','Color',plotcolsummary,'LineWidth',2);
    
    Xlim=xlim;
    Ylim=ylim;
    
    text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-(Ylim(2)-Ylim(1))/10,['p=' num2str(SummaryStats.p)],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
    text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-2*(Ylim(2)-Ylim(1))/10,['R2=' num2str(SummaryStats.R2)],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
    text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-3*(Ylim(2)-Ylim(1))/10,['slope=' num2str(SummaryStats.slope) '(' num2str(SummaryStats.slopeCI) ')'],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
    text(Xlim(2)-(Xlim(2)-Xlim(1))/2,Ylim(2)-4*(Ylim(2)-Ylim(1))/10,['int=' num2str(SummaryStats.intercept) '(' num2str(SummaryStats.interceptCI) ')'],'FontSize',13,'FontWeight','bold','Color',plotcolsummary);
end


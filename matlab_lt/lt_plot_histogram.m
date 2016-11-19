function [Ybinned, Xcenters, hbar]=lt_plot_histogram(Y, Xcenters, plotON, pdfON, FaceAlpha, StairsON, color)
%% LT 8/13/15 - enter data - makes histogram and plots bar plot
% Y = data
% Xcenters = bin centers OPTIONAL (if don't give, will choose based on
% range of data)
% FaceAlpha = 0.3 transparency
% StairsON=1,. plots stairs. otherwise plots bars

% outputs:
% Ybinned = binned data
% hbar = handle of barplot


if ~exist('Xcenters', 'var');
     minbin=min(Y)-std(Y)/3;
    maxbin=max(Y)+std(Y)/3;
    
    Nbins=ceil(numel(Y)/4);
    if Nbins<8;
        Nbins=8;
    end
    Xcenters=linspace(minbin, maxbin, Nbins);
    
elseif isempty(Xcenters);
    minbin=min(Y)-std(Y)/3;
    maxbin=max(Y)+std(Y)/3;
    
    Nbins=ceil(numel(Y)/4);
    if Nbins<8;
        Nbins=8;
    end
    Xcenters=linspace(minbin, maxbin, Nbins);
end

if ~exist('color','var');
    color='k';
end

if ~exist('plotON', 'var');
    plotON=1;
end

if ~exist('pdfON', 'var');
    pdfON=0;
end

if ~exist('FaceAlpha', 'var');
    FaceAlpha=0.3;
end
if isempty(FaceAlpha)
    FaceAlpha=0.3;
end

if ~exist('StairsON', 'var');
    StairsON=0;
end


% get binned data
[Ybinned, Xcenters]=hist(Y, Xcenters);

% get pdf 
if pdfON==1;
    Ybinned=Ybinned./sum(Ybinned);
end


% Plot
if plotON==1;
    if StairsON==0
    hbar=lt_plot_bar(Xcenters, Ybinned, {'FaceAlpha',FaceAlpha, 'Color',color});
    
    else
       binwidth=Xcenters(2)-Xcenters(1); % go from center to left edge
%        inds=find(Ybinned>0);  % remove x with no data, or else will get line at 0
        hbar=stairs(Xcenters-binwidth/2, Ybinned);
        set(hbar, 'LineWidth' , 3, 'Color',color)
%         Ylim=ylim; % move slightly up so line at zero is not visible
%         ylim([0.01 Ylim(2)])
    end
        
    if pdfON==1;
        ylabel('prob density');
    else
        ylabel('count');
    end
    
    % plot mean + sem
    Ymean=mean(Y);
    Ysem=lt_sem(Y);

    Ylim=ylim;
    Yrange=Ylim(2)-Ylim(1);
    plot(Ymean, Ylim(2)-Yrange/10, 'd', 'MarkerSize',8, 'Color',color, 'MarkerFaceColor',color);
    line([Ymean-Ysem, Ymean+Ysem], [Ylim(2)-Yrange/10 Ylim(2)-Yrange/10], 'Color', color, 'LineWidth',2);

    
else
    hbar=[];
end


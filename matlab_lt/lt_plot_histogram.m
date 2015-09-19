function [Ybinned, Xcenters, hbar]=lt_plot_histogram(Y, Xcenters, plotON, pdfON)
%% LT 8/13/15 - enter data - makes histogram and plots bar plot
% Y = data
% Xcenters = bin centers OPTIONAL (if don't give, will choose based on
% range of data)

% outputs:
% Ybinned = binned data
% hbar = handle of barplot


if ~exist('Xcenters', 'var');
     minbin=min(Y)-std(Y)/3;
    maxbin=max(Y)+std(Y)/3;
    
    Nbins=ceil(numel(Y)/4);
    if Nbins<2;
        Nbins=2;
    end
    Xcenters=linspace(minbin, maxbin, Nbins);
elseif isempty(Xcenters);
         minbin=min(Y)-std(Y)/3;
    maxbin=max(Y)+std(Y)/3;
    
    Nbins=ceil(numel(Y)/4);
    if Nbins<2;
        Nbins=2;
    end
    
    Xcenters=linspace(minbin, maxbin, Nbins);
   
end

if ~exist('plotON', 'var');
    plotON=1;
end

if ~exist('pdfON', 'var');
    pdfON=0;
end


% get binned data
[Ybinned, Xcenters]=hist(Y, Xcenters);

% get pdf 
if pdfON==1;
    Ybinned=Ybinned./sum(Ybinned);
end


% Plot
if plotON==1;
hbar=lt_plot_bar(Xcenters, Ybinned);

if pdfON==1;
    ylabel('prob density');
else
    ylabel('count');
end
end


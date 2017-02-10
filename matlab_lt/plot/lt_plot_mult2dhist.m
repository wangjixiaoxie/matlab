function lt_plot_mult2dhist(Y, x, Xcenters, colorplot, transparentOn)
%NOTE: colorplot must be vector

% Y = []; % 2d matrix, each col will be a different histogram, plotted in order
% x = []; % array, with one val for each column

% xaxis = variable of interest
% zaxis = probability
% yaxis = whatever you put in "x";

%%

% 1) -- determine appropriate xbin edges
Ytmp = reshape(Y, numel(Y), 1); % convert to vector

if ~exist('Xcenters', 'var');
[~, Xcenters] = lt_plot_histogram(Ytmp, '', 0, 0, 0, 0, '');
tmp = ceil(size(Y,2)/2);
Xcenters = Xcenters(1:tmp:end);
end

if ~exist('colorplot', 'var')
    colorplot = [1 0 0];
end

if ~exist('transparentOn', 'var')
    transparentOn = 0;
end
   


%%
% 1) collect vertices for each column
Xvert = {};
Yvert = {};
Zvert = {};
for i=1:size(Y,2)
    hfig  = figure; hold on;
    y = Y(:,i);
    
    if (0)
        hist(y, Xcenters);
        
        h = get(gca, 'Children');
        tmp = get(h, 'Vertices');
        
        Xvert{i} = tmp(:,1);
        Zvert{i} = tmp(:,2);
        Yvert{i} = x(i)*ones(size(tmp, 1));
    else
        [~,~,hbar] = lt_plot_histogram(y, Xcenters, 1, 1, 1, 1, 'r');
        h = get(hbar, 'Children');
        Xvert{i} = get(h, 'XData');
        Xvert{i} = [Xvert{i}(1) Xvert{i} Xvert{i}(end)]; % to make sure ends at z = 0
        Zvert{i} = get(h, 'YData');
        Zvert{i} = [0 Zvert{i} 0];
        Yvert{i} = x(i)*ones(1, length(Xvert{i}));

%         Xvert{i} = get(h, 'XData');
%         Zvert{i} = get(h, 'YData');
%         Yvert{i} = x(i)*ones(1, length(Xvert{i}));

    end
    
    close(hfig);
end


% lt_figure; hold on;
plotcols = lt_make_plot_colors(length(Xvert), 1, colorplot);
for i=1:length(Xvert)
    % p = patch(Xvert{i}, Yvert{i}, Zvert{i}, 'FaceColor', 'k');
    p = patch(Xvert{i}', Yvert{i}', Zvert{i}', 'k');
    if transparentOn == 0 
    set(p, 'FaceColor', plotcols{i}, 'EdgeColor', 'k');
    else
    set(p, 'FaceColor', 'none', 'EdgeColor', plotcols{i}, 'LineWidth', 2);
    end        
    % set(p, 'FaceColor', 'none', 'EdgeColor', 'k');
end

view(3)

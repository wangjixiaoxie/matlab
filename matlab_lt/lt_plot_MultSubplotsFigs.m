%% LT 7/20/15 - gives you correct figure and subplot number, so you can put multiple subplots over multiple figs

function [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

% To run:

% Params:
% figcount = 1 % Put this OUTSIDE the for loop
% SubplotsPerFig - how many subplots per fig?
% subplotrows and subplotcols should match SubplotsPerFig
% fignums_alreadyused=[];
% hfigs=[];

% EXAMPLE:
% figcount=1;
% subplotrows=2;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];

% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);


% put this function where you would normally put subplot
% put the params above before the entire for loop.

% This code automatically pulls up the correct figure and subplot.  (i.e.
% no need to call figure(fignum) and subplot(subplotnum);

%% PARAMS

% fignums_alreadyused=[];
% hfigs=[];
SubplotsPerFig=subplotrows*subplotcols;

%% RUN

fignum=ceil(count/SubplotsPerFig);

% Keep track of all fignums alraedy used
fignums_alreadyused=[fignums_alreadyused fignum];

% if this is new figure, get a new handle
if length(fignums_alreadyused)==1 || fignums_alreadyused(end)~=fignums_alreadyused(end-1);
    % then is new figure
    hfigs(fignum)=figure; hold on;
    lt_plot_format;
else
    % is old figure, just reopen it
    figure(hfigs(fignum));
end

% Which subplot is this?
subplot_num=count-(fignum-1)*SubplotsPerFig;
hsplot=lt_subplot(subplotrows, subplotcols, subplot_num); hold on;
set(gca, 'YColor', [0.5 0.5 0.5])
set(gca, 'XColor', [0.5 0.5 0.5])
subplotsqueeze(gca, 1.2);
% Update count
count=count+1;





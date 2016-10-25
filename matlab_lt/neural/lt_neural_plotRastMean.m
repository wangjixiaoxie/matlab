%%  lt 8/19/16 - plots raster and (pltos and outputs) mean +/- SE (overlapping)

function [xbin, ymean, ysem, ymean_hz, ysem_hz, ystd, ystd_hz] = lt_neural_plotRastMean(Yspks, window, windshift, plotRastOn, plotcol)


% %% ==== [for converting structure to cell array]
%
% % - convert to cell format
% numsegs=length(SegmentsExtract);
% k=1;
% for i=1:numsegs
%
%         inds=find([SegmentsExtract(i).spk_Clust]==k);
%
%     spktimes=SegmentsExtract(i).spk_Times(inds);
%
%     Yspks{i}=spktimes;
% end

%% === input params
% plotcol='k';
% Yspks; a cell, where each ind is a vector of spike times. each ind is a
% different trial.
%
% window=0.015; sec for window.
% windshift=0.002; % how much to slide window
% plotRastOn=1; % then plots rasters

% out
% ymean vector of mean spks/bin
% ysem vector of se for each bin
% xbin= time of center of each bin
%% === PLOT RASTERS ON ACTIVE AXES
numsegs=length(Yspks); % num trials

if plotRastOn==1
    
    for i=1:numsegs
        
        if isempty(Yspks{i})
            continue
        else
            for ii=1:length(Yspks{i})
                spktime=Yspks{i}(ii);
                line([spktime spktime], [i-0.4 i+0.4], 'Color', plotcol, 'LineWidth',0.5);
            end
        end
    end
end

%% === OUTPUT RUNNING MEAN

xmax=max(cellfun(@max, Yspks(~cellfun(@isempty, Yspks))));

xstarts=0:windshift:xmax-window;
xends=window:windshift:xmax;

Xall=[]; % trial x bins
Yall=[];
for i=1:numsegs
    
    spktimes=Yspks{i};
    
    Xseg=[];
    Yseg=[];
    for j=1:length(xstarts)
        % for this window, get count
        wmin=xstarts(j);
        wmax=xends(j);
        
        y=sum(spktimes>wmin & spktimes<wmax);
        
        Xseg=[Xseg mean([wmin wmax])]; % x val is bincenter
        Yseg=[Yseg y];
    end
    
    Xall=[Xall; Xseg];
    Yall=[Yall; Yseg];
end

% -- get mean and sem across trials
ymean=mean(Yall, 1);
ystd=std(Yall, 0, 1);
ysem=lt_sem(Yall);
xbin=Xall(1,:);

% - convert to spike rate (assumes that binwidth is in sec)
ymean_hz=ymean./window;
ysem_hz=ysem./window;
ystd_hz=ystd./window;




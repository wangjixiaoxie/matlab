function lt_neural_v2_CTXT_PlotGeneral_M2(CLASSEScompiled, plotstat)
%% lt 8/15/17 - takes output from lt_neural_v2_CTXT_PlotGeneral_M2 and plots

%% ---

numanalys = length(CLASSEScompiled.analynum); % actually is frwindow sizes

%%
% ---- collect params across all experiments
tmp.AllFRwindsize = [];
tmp.AllFRtwindow = [];
tmp.FRbinsize = [];
tmp.algnonset = [];
tmp.algnsyl = [];


for i = 1:numanalys % fr window size
    
    numiterations = length(CLASSEScompiled.analynum(i).iterationnum);
    frwindowsize = diff(CLASSEScompiled.analynum(i).iterationnum(1).params.ClassGeneral.frtimewindow);
    
    for ii=1:numiterations % combination of fr bin size and window location
            
        frbinsize = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frbinsize;
        frtimewindow = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frtimewindow;

        algnsyl = CLASSEScompiled.analynum(i).iterationnum(ii).params.alignWhichSyl;
        algnonset = CLASSEScompiled.analynum(i).iterationnum(ii).params.alignOnset;
        

tmp.AllFRwindsize = [tmp.AllFRwindsize; frwindowsize];
tmp.AllFRtwindow = [tmp.AllFRtwindow; frtimewindow];
tmp.FRbinsize = [tmp.FRbinsize; frbinsize];
tmp.algnonset = [tmp.algnonset; algnonset];
tmp.algnsyl = [tmp.algnsyl; algnsyl];
        
    end
end

assert(length(unique(tmp.algnsyl)) ==1, 'problem, diff syl algn experiments included ...');
assert(length(unique(tmp.algnonset)) ==1, 'problem, diff algn onset experiments included ...');

disp(['FR window sizes used: ' num2str(unique(tmp.AllFRwindsize)')])
disp(['FR binsizes used: ' num2str(unique(tmp.FRbinsize)')])

SylAlignNum = unique(tmp.algnsyl);
SylAlignOnset = unique(tmp.algnonset);

clear tmp;

%% PLOT
AllFRwindowsize = [];
AllFRwindow_middle = [];
AllFRbinsize = [];
AllPerformance = [];
AllNumctxts = [];

% NOTE: below collects all syl contours across everything (including all
% branches and numctxts) - should change to at least split by branch type
% (or at least num ctxts);
        AllSylContours = [];
        AllSylContours_xtimes = [];
        
for i = 1:numanalys % fr window size
    
    numiterations = length(CLASSEScompiled.analynum(i).iterationnum);
    frwindowsize = diff(CLASSEScompiled.analynum(i).iterationnum(1).params.ClassGeneral.frtimewindow);
    
    for ii=1:numiterations % combination of fr bin size and window location
        
        frbinsize = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frbinsize;
        frtimewindow = CLASSEScompiled.analynum(i).iterationnum(ii).params.ClassGeneral.frtimewindow;
        windowmid = mean(frtimewindow);
        
        accuracy = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllAccuracy];
        sensitivity = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllSensitivity];
        piyactual = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllPiYActualMean];
        
        Numctxts = [CLASSEScompiled.analynum(i).iterationnum(ii).allbranches.AllNumCtxts];
        
%         xbins = frtimewindow(1):frbinsize:frtimewindow(end);
%         xbins = mean([xbins(1:end-1)' xbins(2:end)'], 2);
% 
        
        % === output
        if strcmp(plotstat, 'accuracy')
            AllPerformance = [AllPerformance accuracy];
        elseif strcmp(plotstat, 'sensitivity')
            AllPerformance = [AllPerformance sensitivity];
        elseif strcmp(plotstat, 'pihatmean')
            AllPerformance = [AllPerformance piyactual];
        else
            disp('problem!!!')
        end
        
        AllNumctxts = [AllNumctxts Numctxts];
        
        % -- expand to sample size
        AllFRbinsize = [AllFRbinsize ones(size(accuracy))*frbinsize];
        AllFRwindow_middle = [AllFRwindow_middle ones(size(accuracy))*windowmid];
        AllFRwindowsize = [AllFRwindowsize ones(size(accuracy))*frwindowsize];
        
        
        % --- save mean syl traces
        % for now, just take mean over all branches
        numbranches = length(CLASSEScompiled.analynum(i).iterationnum(ii).allbranches);
        motifpredur = CLASSEScompiled.analynum(i).iterationnum(ii).params.motifpredur;
        
        SylContourstmp = [];
           
        for bb=1:numbranches
            
            sylcontours = CLASSEScompiled.analynum(i).iterationnum(ii).allbranches(bb).AllBranchSylContours;
          
            SylContourstmp = [SylContourstmp; mean(sylcontours,1)];

        end
        
        AllSylContours = [AllSylContours mean(SylContourstmp,1)];
        xtimes = (1:size(sylcontours,2))/1000 - motifpredur;
        AllSylContours_xtimes = [AllSylContours_xtimes xtimes];

        
    end
end

[SylContourMean, ystd] = grpstats(AllSylContours, AllSylContours_xtimes, {'mean', 'std'});
SylContourMean_times = unique(AllSylContours_xtimes);

%% =================== PLOT
% Separate plots for each: 0) value of numctxts 1) binsize, 2) Fr window size

figcount=1;
subplotrows=7;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

binsizelist = unique(AllFRbinsize);
numctxtslist = unique(AllNumctxts);
frwindowlist = unique(AllFRwindowsize);
hsplots = [];


for fw = frwindowlist
    for k=numctxtslist
        
        for bsize = binsizelist;
            
            inds = find(AllFRbinsize==bsize & AllFRwindowsize==fw & AllNumctxts==k);
            
            if isempty(inds)
                continue
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['wsize:' num2str(fw) ', bsize:' num2str(bsize) ', Nctxt:' num2str(k)]);
            
            xtimes = AllFRwindow_middle(inds);
            yvals = AllPerformance(inds);
            
            plot(xtimes, yvals, 'x', 'Color', [0.6 0.6 0.6]);
            
            [ymean, ysem] = grpstats(yvals, xtimes, {'mean', 'sem'});
            xtimes_grp = unique(xtimes);
            lt_plot(xtimes_grp, ymean, {'Errors', ysem, 'Color', 'k'});
            axis tight
            
            % -- line indicating null value for this number of contexts
            line(xlim, [1/k 1/k], 'Color', 'b');
            
            % --
            plot(SylContourMean_times, SylContourMean, '-r');
        end
    end
end
linkaxes(hsplots, 'xy');


%%
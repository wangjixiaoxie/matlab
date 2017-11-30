function lt_neural_v2_CTXT_BRANCH_PlotByBranchID(ALLBRANCH)
%% lt 11/29/17 - plots, sepaated by unique branch points (i.e. regexpstr)


%% ===== organize data by unique branches

DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH);

%% ===== FOR EACH BIRD PLOT EACH BRANCH, SEPARATING BY BRAIN REGION
numbirds = length(DATSTRUCT_BYBRANCH.bird);

locationstoplot = {'LMAN', 'RA'};


for i=1:numbirds
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        for loc = locationstoplot
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            
            % -------------- PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([thisbranch{1} ' [' loc{1} ']']);
                        
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            
            plot(xdecode, ydecode, '-', 'Color', [0.4 0.4 0.4]);
            
            % -- plot mean
            if size(ydecode,1)>1
            ymean = mean(ydecode,1);
            ysem = lt_sem(ydecode);
            shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
            end
            
            % -- lines
            line([0 0], ylim);
        end
        
    end
    axis tight
    linkaxes(hsplots, 'xy');
end
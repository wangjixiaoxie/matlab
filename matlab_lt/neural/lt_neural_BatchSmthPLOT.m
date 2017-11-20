function lt_neural_BatchSmthPLOT(DATSTRUCT, chanstoplot, motifstoplot, ...
    sylstoalign, pretime, posttime, plotRaw, fs)

DirFieldList = fieldnames(DATSTRUCT);
assert(length(DirFieldList)==2, 'problem - code not ready for >2 fileds yet');

FRwind = [-0.035 0.02];


%%
% ============================ PLOT RESULTS
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% ==== three subplots for each motif, dir, undir, and overlay dir and undir
nummotifs = length(DATSTRUCT.(DirFieldList{1}).motifnum);
numchannels = length(chanstoplot);


% ====== to collect correlation (of mean sm rate across trials)
for i=1:nummotifs
    motifname = DATSTRUCT.(DirFieldList{1}).motifnum(i).motifname;
    
    %     = ~cellfun('isempty', DATSTRUCT_DIR.motifnum(i).DatAll)
    
    for cc=chanstoplot
        hsplots = [];
        % ----- UNDIR
        Ymean_UNDIR = [];
        Ysem_UNDIR = [];
        dirfield = DirFieldList{1};
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['[' dirfield ']' motifname '-ch' num2str(cc)]);
        hsplots = [hsplots hsplot];
        datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{cc};
        t = DATSTRUCT.(dirfield).motifnum(i).t;
        if ~isempty(datmat)
        % -- plot raw
        plot(t, datmat, 'Color', [0.7 0.7 0.7]);
        % -- overlay median, with 75th tiles
        datmedian = median(datmat,1);
        datCI = prctile(datmat, [75 25]);
        datSEM = lt_sem(datmat);
        shadedErrorBar(t, datmedian, datSEM, {'Color','r'}, 1);
        line([pretime pretime],ylim);
        lt_plot_zeroline;
        axis tight;
        
        Ymean_UNDIR = datmedian;
        Ysem_UNDIR = datSEM;
        end
        
        % ---- DIR
        Ymean_DIR = [];
        Ysem_DIR = [];
        dirfield = DirFieldList{2};
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[' dirfield ']' motifname '-ch' num2str(cc)]);
        datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{cc};
        t = DATSTRUCT.(dirfield).motifnum(i).t;
        % -- plot raw\
        if ~isempty(datmat)
        plot(t, datmat, 'Color', [0.7 0.7 0.7]);
        % -- overlay median, with 75th tiles
        datmedian = median(datmat,1);
        datCI = prctile(datmat, [75 25]);
        datSEM = lt_sem(datmat);
        shadedErrorBar(t, datmedian, datSEM, {'Color','r'}, 1);
        line([pretime pretime],ylim);
        lt_plot_zeroline;
                axis tight;

        Ymean_DIR = datmedian;
        Ysem_DIR = datSEM;
        end
        
        % ----- COMBINED
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[DIR(r) AND UNDIR(k)]' motifname '-ch' num2str(cc)]);
        
        if ~isempty(Ymean_UNDIR)
        shadedErrorBar(t, Ymean_UNDIR, Ysem_UNDIR, {'Color','k'}, 1);
        end
        if ~isempty(Ymean_DIR)
        shadedErrorBar(t, Ymean_DIR, Ysem_DIR, {'Color','r'}, 1);
        line([pretime pretime],ylim);
        lt_plot_zeroline;
                axis tight;
        end

        % ---
        linkaxes(hsplots, 'xy');
        
    end
end


%% ================= correlations across trials
% collect all pairwise correlations between channels

figcount=1;
subplotrows=3;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

windwind = pretime+FRwind; % convert to time rel onset

            % ====== to collect correlation (of mean sm rate across trials)
for i=1:nummotifs
    motifname = DATSTRUCT.(DirFieldList{1}).motifnum(i).motifname;
    
    %     = ~cellfun('isempty', DATSTRUCT_DIR.motifnum(i).DatAll)
    
    for j =1 :length(chanstoplot)
        
        for jj=j+1:length(chanstoplot)
            
            chan1 = chanstoplot(j);
            chan2 = chanstoplot(jj);
            
            % ========== UNDIR
            dirfield = DirFieldList{1};
            plotcol = 'k';
            
            FRtrials_1 = [];
            % -- 1st chan 1
%             if any(chan1 == [9 14]); % then in is LMAN
%                 windwind = [-0.035 0.005];
%             elseif chan1 == 21
%                 windwind = [-0.015 0.025];
%             end
            t = DATSTRUCT.(dirfield).motifnum(i).t;
            inds = t>=windwind(1) & t<windwind(2);
            datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{chan1};
            if ~isempty(datmat)
            FRtrials_1 = mean(datmat(:,inds), 2); % collect one mean FR for each trial
%             FRtrials_1 = FRtrials_1 - mean(FRtrials_1); % subtract mean
            FRtrials_1 = (FRtrials_1 - mean(FRtrials_1))./std(FRtrials_1); % subtract mean
            else
                FRtrials_1 = [];
            end
            
            FRtrials_2 = [];
            % -- 2nd chan
            t = DATSTRUCT.(dirfield).motifnum(i).t;
            inds = t>=windwind(1) & t<windwind(2);
            datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{chan2};
            if ~isempty(datmat)
            FRtrials_2 = mean(datmat(:,inds), 2);
            FRtrials_2 = (FRtrials_2 - mean(FRtrials_2))./std(FRtrials_2); % subtract mean
            else
                FRtrials_2 = [];
            end
            
           % ======== plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['[' dirfield ']-' motifname '-ch' num2str(chan1) ' vs ' num2str(chan2)]);
            xlabel(['ch' num2str(chan1)]);
            ylabel(['ch' num2str(chan2)]);
            lt_regress(FRtrials_2, FRtrials_1, 1, 0, 1, 1, plotcol);
            
            
            
            % ========== DIR
            dirfield = DirFieldList{2};
            plotcol = 'r';
            
            FRtrials_1 = [];
            % -- 1st chan 1
            t = DATSTRUCT.(dirfield).motifnum(i).t;
            inds = t>=windwind(1) & t<windwind(2);
            datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{chan1};
            FRtrials_1 = mean(datmat(:,inds), 2); % collect one mean FR for each trial
            FRtrials_1 = (FRtrials_1 - mean(FRtrials_1))./std(FRtrials_1); % subtract mean
            
            FRtrials_2 = [];
            % -- 2nd chan
            t = DATSTRUCT.(dirfield).motifnum(i).t;
            inds = t>=windwind(1) & t<windwind(2);
            datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{chan2};
            FRtrials_2 = mean(datmat(:,inds), 2);
            FRtrials_2 = (FRtrials_2 - mean(FRtrials_2))./std(FRtrials_2); % subtract mean
            
           % ======== plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(['[' dirfield ']-' motifname '-ch' num2str(chan1) ' vs ' num2str(chan2)]);

            xlabel(['ch' num2str(chan1)]);
            ylabel(['ch' num2str(chan2)]);
            lt_regress(FRtrials_2, FRtrials_1, 1, 0, 1, 1, plotcol);
            
            
        end
                
    end
end

linkaxes(hsplots, 'xy');


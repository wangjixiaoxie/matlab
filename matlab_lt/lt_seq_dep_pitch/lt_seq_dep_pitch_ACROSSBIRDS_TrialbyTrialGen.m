function lt_seq_dep_pitch_ACROSSBIRDS_TrialbyTrialGen(TrialStruct, ParamsTrial, ...
    SeqDepPitch_AcrossBirds, plotSongBySong, plotRaw)
%% to do (11/28/17 )
1) compare raw FF to cross cov functions for each experiment
2) only look at experiments with strong learning
3) only look at early learning + last baseline day (combine into one vector
allows 


%%
% -- for interpolation
xbinsize = 1; % in units minutes
interpmethod = 'pchip';


%% ==== stuff

% convert interpolation binsize to days
xbinsize = xbinsize/(24*60);
NumBirds = length(TrialStruct.birds);

%% =========== get ffvals by song for each syl

for i=1:NumBirds
    
    numexpts = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:numexpts
        
if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp('SKIP - no songs');
            continue
        end
                numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        
        % ---- get time of all unique songs for this experiment
        %         tvalsAll = [];
        %         for j=1:numsyls
        %            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals_datenum;
        %            tvalsAll = [tvalsAll; tvals];
        %         end
        %         tvalsUnique = unique(tvalsAll);
        %
        %         accumuarray(subs, val, [], @mean);
        
        for j=1:numsyls
            
            tvalsdnum = TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals_datenum;
            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals;
            ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).FFvals;
            
            % ------------ convert to song by song FF
            [~, indstmp] = sort(tvalsdnum);
            tvalsdnum = tvalsdnum(indstmp);
            tvals = tvals(indstmp);
            ffvals = ffvals(indstmp);
            
            ffvals = grpstats(ffvals, tvals);
            tvals = unique(tvals);
            tvalsdnum = unique(tvalsdnum);
            assert(length(tvals) == length(tvalsdnum), 'asdfas');
            
            % ======================= SAVE
            TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_Tvals_datenum = tvalsdnum;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_FFvals = ffvals;
            TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_Tvals = tvals;
        end
        
    end
end



%% ========= 1) EACH SONG AS ONE DATAPOINT - SONG BY SONG GENERALIZATION
if plotRaw==1
for i=1:NumBirds
    
    numexpts = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:numexpts
        
        
if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp('SKIP - no songs');
            continue
end

figcount=1;
        subplotrows=5;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        
        % ------- day parameters
        baseDays = TrialStruct.birds(i).exptnum(ii).BaseDays;
        
        
        
        % ================ first plot targ
        targind = [TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget];
        targind = logical(targind);
        if plotSongBySong==1
            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(targind).songbysong_Tvals;
            ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(targind).songbysong_FFvals;
        else
            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(targind).Tvals;
            ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(targind).FFvals;
        end
        syl = TrialStruct.birds(i).exptnum(ii).sylnum(targind).syl;
        
        % --- subtract baseline
        ffvals = ffvals - mean(ffvals(tvals<max(baseDays)+1));
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots =[hsplots hsplot];
        title(syl);
        plot(tvals, ffvals, 'ok');
        
        % --- lines
        line([baseDays(end)+1 baseDays(end)+1], ylim);
        lt_plot_zeroline;
        
        % ==================== plot other syls
        for j=1:numsyls
            
            % --- skip if is targ
            if TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_istarget==1
                continue
            end
            syl = TrialStruct.birds(i).exptnum(ii).sylnum(j).syl;
            
            if plotSongBySong==1
                tvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_Tvals;
                ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_FFvals;
            else
                tvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals;
                ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).FFvals;
            end
            
            issimilar = TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_similar;
            if issimilar==1
                plotcol = 'b';
            else
                plotcol = 'r';
            end
            
            % --- subtract baseline
            ffvals = ffvals - mean(ffvals(tvals<max(baseDays)+1));
            
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(syl);
            hsplots =[hsplots hsplot];
            plot(tvals, ffvals, 'o', 'Color', plotcol);
            line([baseDays(end)+1 baseDays(end)+1], ylim);
            lt_plot_zeroline;
            
        end
        linkaxes(hsplots, 'xy');
    end
end
end

%% =================== correlation of between syls, using song by song dat
AllXcov = struct;
AllXcov.sametype.base = [];
AllXcov.sametype.WNon = [];
AllXcov.diff.base = [];
AllXcov.diff.WNon = [];

plotOn = 0;
            windxcov = 10;
            binxcov = 1;
subtractrunningmean = 1;

for i=1:NumBirds
    
    numexpts = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:numexpts
        
        
        if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp('SKIP - no songs');
            continue
        end
        
        if plotOn==1
                figcount=1;
                subplotrows=5;
                subplotcols=2;
                fignums_alreadyused=[];
                hfigs=[];
                hsplots = [];
        end
        
        numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        
        % ------- day parameters
        baseDays = TrialStruct.birds(i).exptnum(ii).BaseDays;
        
        
        % ================ first get targ
        targind = [TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget];
        targind = logical(targind);
%         syl = TrialStruct.birds(i).exptnum(ii).sylnum(targind).syl;
                
        tvalsTarg = TrialStruct.birds(i).exptnum(ii).sylnum(targind).songbysong_Tvals;
        ffvalsTarg = TrialStruct.birds(i).exptnum(ii).sylnum(targind).songbysong_FFvals;
        
        % ------
        if subtractrunningmean==1
            binsize = 15; % must be odd
            tmp = lt_running_stats(ffvalsTarg, binsize);
            ffmean = [ones(1,(binsize-1)/2)*tmp.Mean(1) tmp.Mean ones(1,(binsize-1)/2)*tmp.Mean(end)];
            ffvalsTarg = ffvalsTarg - ffmean';
        end
        
        % ==================== compare other syls to targ
        for j=1:numsyls
            
            if TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_istarget==1
                continue
            end
            
            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_Tvals;
            ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).songbysong_FFvals;
            syl = TrialStruct.birds(i).exptnum(ii).sylnum(j).syl;
            issimilar = TrialStruct.birds(i).exptnum(ii).sylnum(j).INFO_similar;
            
    % ----
            if subtractrunningmean==1
            binsize = 15; % must be odd
            tmp = lt_running_stats(ffvals, binsize);
            ffmean = [ones(1,(binsize-1)/2)*tmp.Mean(1) tmp.Mean ones(1,(binsize-1)/2)*tmp.Mean(end)];
            ffvals = ffvals - ffmean';
        end

            
            % ------------------------- pull out songs for which they both
            % have data
            [~, ind1, ind2] = intersect(tvalsTarg, tvals);
            
            Ttarg = tvalsTarg(ind1);
            Ftarg = ffvalsTarg(ind1);
            
            Tsyl = tvals(ind2);
            Fsyl = ffvals(ind2);
            assert(length(ind1)==length(ind2), 'asdfasd');
            
            
            % ========================= correlate with target
            if plotOn==1
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            if issimilar==1
                title([syl], 'Color', 'b');
                typefname = 'sametype';
            else
                title([syl], 'Color', 'k');
                typefname = 'diff';
            end
            xlabel(['targ <--> ' syl]);
            end
            
            % ----- for output
            if issimilar==1
                typefname = 'sametype';
            else
                typefname = 'diff';
            end
            
            % ############ BASELINE
            f1 = Ftarg(Ttarg<baseDays(end)+1); 
            f2 = Fsyl(Tsyl<baseDays(end)+1);
            
            [cc, lag] = xcov(f1, f2, windxcov/binxcov, 'coeff');

            if plotOn==1
            plot(lag*binxcov, cc, '-ok');
            lt_plot_zeroline;
            end
            
            % --- colect
            AllXcov.(typefname).base = [AllXcov.(typefname).base; cc'];
            
            
            % ############ TRAINING
            f1 = Ftarg(Ttarg>baseDays(end)+1); 
            f2 = Fsyl(Tsyl>baseDays(end)+1);
            
            [cc, lag] = xcov(f1, f2, windxcov/binxcov, 'coeff');

            if plotOn==1
            plot(lag*binxcov, cc, '-or');
            lt_plot_zeroline;
            end
            
            % --- colect
            AllXcov.(typefname).WNon = [AllXcov.(typefname).WNon; cc'];
            
           
        end
    end
end


% ======================== PLOT SUMMARY
lt_figure; hold on;

% %%%%%%%%%%%%%%%%%%%%%% same
lt_subplot(1,2,1); hold on;
title('same');
% --- base
plotcol = 'k';
ccmean = mean(AllXcov.sametype.base,1);
ccsem = lt_sem(AllXcov.sametype.base);
shadedErrorBar(lag, ccmean, ccsem, {'Color', plotcol}, 1);

% --- wn on
plotcol = 'r';
ccmean = mean(AllXcov.sametype.WNon,1);
ccsem = lt_sem(AllXcov.sametype.WNon);
shadedErrorBar(lag, ccmean, ccsem, {'Color', plotcol}, 1);

lt_plot_zeroline;
ylim([-0.2 0.8]);

% %%%%%%%%%%%%%%%%%%%%%%%%% diff
lt_subplot(1,2,2); hold on;
title('diff');
% --- base
plotcol = 'k';
ccmean = mean(AllXcov.diff.base,1);
ccsem = lt_sem(AllXcov.diff.base);
shadedErrorBar(lag, ccmean, ccsem, {'Color', plotcol}, 1);

% --- wn on
plotcol = 'r';
ccmean = mean(AllXcov.diff.WNon,1);
ccsem = lt_sem(AllXcov.diff.WNon);
shadedErrorBar(lag, ccmean, ccsem, {'Color', plotcol}, 1);
lt_plot_zeroline;
ylim([-0.2 0.8]);

%% ================== 1) interpolate to get continuous signal (can't use trials since not fully labeled).


for i=1:NumBirds
    
    numexpts = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:numexpts
        
        
        numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        
        % === for each syl, extract timecourse
        % -- get x base
        func = @(X) max(X);
        maxtime = cellfun(func, {TrialStruct.birds(i).exptnum(ii).sylnum.Tvals});
        func = @(X) min(X);
        mintime = cellfun(func,  {TrialStruct.birds(i).exptnum(ii).sylnum.Tvals});
        
        xvals = mintime-xbinsize:xbinsize:maxtime+xbinsize;
        
        YvalsAll = nan(numsyls, length(xvals)); % syls x timepts
        
        for j=1:numsyls
            
            tvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).Tvals;
            ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(j).FFvals;
            
            % ------------ convert to song by song FF
            [~, indstmp] = sort(tvals);
            tvals = tvals(indstmp);
            ffvals = ffvals(indstmp);
            
            ffvals = grpstats(ffvals, tvals);
            tvals = unique(tvals);
            
            %         lt_figure; hold on;
            %
            %         subplot(411); title('linear');
            %         yvals = interp1(tvals, ffvals, xvals, 'linear');
            %         plot(xvals, yvals, '-r', tvals, ffvals, 'ko');
            %
            %                 subplot(412); title('nearest');
            %         yvals = interp1(tvals, ffvals, xvals, 'nearest');
            %         plot(xvals, yvals, '-r', tvals, ffvals, 'ko');
            %
            %
            %         subplot(413); title('cubic');
            yvals = interp1(tvals, ffvals, xvals, 'pchip');
            %         plot(xvals, yvals, '-r', tvals, ffvals, 'ko');
            %
            %         subplot(414); title('spline');
            %         yvals = interp1(tvals, ffvals, xvals, 'spline');
            %         plot(xvals, yvals, '-r', tvals, ffvals, 'ko');
            %         pause
            %         close all
            
            % ================= COLLECT ALL SYLS
            YvalsAll(j,:) = yvals;
            
        end
        assert(~any(isnan(YvalsAll(:))), 'asfasd');
        YvalsAll = int16(YvalsAll);
        
        % ========== STORE
        TrialStruct.birds(i).exptnum(ii).FFinterpAll_sylxtime = YvalsAll;
        TrialStruct.birds(i).exptnum(ii).FFinterpAll_tvals = tvals;
        
    end
end

%% ============= COHERENCY OR CORRELATION AT DIFFERTNT TIME LAGS (BASE VS. TRAINING)

NumBirds = length(TrialStruct.birds);

for i=1:NumBirds
    
    numexpts = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:numexpts
        
        numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        targsyl = TrialStruct.birds(i).exptnum(ii).targsyl;
        targsylind = find(strcmp(TrialStruct.birds(i).exptnum(ii).SylsUnique, targsyl));
        
        YvalsAll = TrialStruct.birds(i).exptnum(ii).FFinterpAll_sylxtime;
        
        yvals_targ = YvalsAll(targsylind, :);
        
        
        for j=1:numsyls
            % for each syl, get
            
        end
        
        
    end
end

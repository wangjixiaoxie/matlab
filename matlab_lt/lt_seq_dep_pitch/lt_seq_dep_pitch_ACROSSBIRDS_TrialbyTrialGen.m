function lt_seq_dep_pitch_ACROSSBIRDS_TrialbyTrialGen(TrialStruct, ParamsTrial, SeqDepPitch_AcrossBirds)

% -- for interpolation
xbinsize = 1; % in units minutes
interpmethod = 'pchip';


%% ==== stuff

% convert interpolation binsize to days
xbinsize = xbinsize/(24*60);



%% ================== 1) interpolate to get continuous signal (can't use trials since not fully labeled). 

NumBirds = length(TrialStruct.birds);

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

function lt_batchsong_plotFF(DATSTRUCT, MotifsToExtract, TrainON, SwitchTimes, subtractMean);

%% LT 11/20/17 - plots

%% INPUTS

% TrainON = '05Nov2017-1125'; WN onset

% MotifsToExtract = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)',  'jbh(h)', '(g)'};


% SwitchTimes = {'05Nov2017-1235', '05Nov2017-1355', '05Nov2017-1548', ...
%     '05Nov2017-1811'}; % will places lines in plot at these times

% subtractMean = 0; (baseline mean)



%%

TrainON_dnum = datenum(TrainON, 'ddmmmyyyy-HHMM');

tval_min = [];
for i=1:length(DATSTRUCT.motif)
   tval_min = min([tval_min min([DATSTRUCT.motif(i).rendnum.datenum_song_SecRes])]);
end

firstday = datestr(tval_min, 'ddmmmyyyy');

%%
figcount=1;
subplotrows=4;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i = 1:length(MotifsToExtract);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(MotifsToExtract{i});
    
    % ####################################### UNDIR
    inds = [DATSTRUCT.motif(i).rendnum.isDIR]==0;
    plotcol = 'k';
    
    % ----------------- RUN
    ffvals = [DATSTRUCT.motif(i).rendnum(inds).ff];
    tvals = [DATSTRUCT.motif(i).rendnum(inds).datenum_song_SecRes];
    
    % ------- subtract baseline FF
    if subtractMean==1
        baseinds = tvals < TrainON_dnum;
        ffvals = ffvals - mean(ffvals(baseinds));
    end
    
    % -- convert tvals to days from start
    tvals = lt_convert_EventTimes_to_RelTimes(firstday, tvals);
    tvals = tvals.FinalValue;
    
    % ----- plot
    plot(tvals, ffvals, 'o', 'Color', plotcol);
    
    % ----- plot day means
    numdays = floor(max(tvals));
    for j=1:numdays
        
        indstmp = floor(tvals)==j;
        tt = tvals(indstmp);
        ff = ffvals(indstmp);
        
        if isempty(ff)
            continue
        end
        
        lt_plot(max(tt)+0.1, mean(ff), {'Errors', lt_sem(ff), 'Color', plotcol});
    end
    
    % ######################### lines
    % --- train onset
    tmp = lt_convert_EventTimes_to_RelTimes(firstday, TrainON_dnum);
    line([tmp.FinalValue tmp.FinalValue], ylim, 'Color', 'r');
    
    for j=1:length(SwitchTimes)
        tmp = lt_convert_EventTimes_to_RelTimes(firstday, datenum(SwitchTimes{j}, 'ddmmmyyyy-HHMM'));
        line([tmp.FinalValue tmp.FinalValue], ylim, 'Color', 'm');
    end
    lt_plot_zeroline;
    
    
    % ####################################### DIR
    inds = [DATSTRUCT.motif(i).rendnum.isDIR]==1;
    if ~any(inds)
        continue
    end
    plotcol = 'b';
    
    % ----------------- RUN
    ffvals = [DATSTRUCT.motif(i).rendnum(inds).ff];
    tvals = [DATSTRUCT.motif(i).rendnum(inds).datenum_song_SecRes];
    
    % ------- subtract baseline FF
    if subtractMean==1
        baseinds = tvals < TrainON_dnum;
        ffvals = ffvals - mean(ffvals(baseinds));
    end
    
    % -- convert tvals to days from start
    tvals = lt_convert_EventTimes_to_RelTimes(firstday, tvals);
    tvals = tvals.FinalValue;
    
    % ----- plot
    lt_plot(tvals, ffvals, {'Color', plotcol});
    %      plot(tvals, ffvals, 'o', 'Color', plotcol);
    
    % ----- plot day means
    numdays = floor(max(tvals));
    for j=1:numdays
        
        indstmp = floor(tvals)==j;
        tt = tvals(indstmp);
        ff = ffvals(indstmp);
        
        if isempty(ff)
            continue
        end
        
        lt_plot(max(tt)+0.15, mean(ff), {'Errors', lt_sem(ff), 'Color', plotcol});
    end
    
    
end
linkaxes(hsplots, 'xy');


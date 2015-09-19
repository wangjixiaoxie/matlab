function [Params, AllDays_RegExpr]=lt_seq_dep_pitch_PlotRegExpr(Params, AllDays_RegExpr ,saveON, LMANon, SuppressQueries)
%% LT 7/13/15 - Starting plots for LMAN data as well.
% LMANon = 1; runs all code only for LMAN data. if not defined, or 0, then
% runs code for non-LMAN data. (runs one or the other, not both)

% SuppressQueries=1 suppresses queries of user - instead will take average
% as input, or something in middle (e.g. for choosing subclass to plot).

%% LT 6/17/15 - Plots data in AllDays_RegExpr


%% PARAMS
if ~exist('SuppressQueries' ,'var');
    SuppressQueries=0;
end


if ~exist('LMANon' ,'var');
    LMANon=0;
end

if ~exist('saveON','var');
    saveON=1;
end

NumDays=length(AllDays_RegExpr.day_data);
NumRegExprClasses=length(Params.RegExpr.expressions);

% --- WN and timeline stuff
FirstDay=Params.SeqFilter.FirstDay;
LastDay=Params.SeqFilter.LastDay;
plotWNdays=Params.PlotRegExpr.plotWNdays;


% Convert WNdays to index days
global WNTimeOnInd % make global, so can use in subfunction below.
global WNTimeOffInd

if Params.PlotRegExpr.plotWNdays==1;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{Params.SeqFilter.WNTimeON});
    WNTimeOnInd=X.JustDays_rel;
    X=lt_convert_EventTimes_to_RelTimes(FirstDay,{Params.SeqFilter.WNTimeOFF});
    WNTimeOffInd=X.JustDays_rel;
end


% Get days to mark, if exist.
global DaysToMarkInds
DaysToMarkInds={};
X={};

if isfield(Params.SeqFilter,'DaysToMark'); % if there are specific days to average over and look at.
    DaysToMark=Params.SeqFilter.DaysToMark;
    for i=1:length(DaysToMark);
        X{i}=lt_convert_EventTimes_to_RelTimes(FirstDay,{DaysToMark{i}});
        DaysToMarkInds{i}=X{i}.JustDays_rel;
    end
end

% UPDATE PARAMS STRUCTURE
Params.PlotRegExpr.WNTimeOnInd=WNTimeOnInd;
Params.PlotRegExpr.WNTimeOffInd=WNTimeOffInd;
Params.PlotRegExpr.DaysToMarkInds=DaysToMarkInds;
Params.PlotRegExpr.LMANon=LMANon;


%% PROCESS DATA - e.g. get relative to baseline, and so forth

% GET SEPARATE STATS FOR EACH CLASS AND SUBCLASS - 
% collect: 1) baseline, 2) day stats, 3) day minus baseline stats

if LMANon==1;   
    data_field='_MUSC';
else
    data_field='';
end


for i=1:NumRegExprClasses;
    
    regexpr_class=Params.RegExpr.expressions{i};
    
    % ---- Get all subclasses
    num_subclasses=length(Params.RegExpr.subexpressions{i});
    
    for ii=1:num_subclasses;
        
        subclass=Params.RegExpr.subexpressions{i}{ii};
        length_subclass=length(subclass);
        
        % =======   BASELINE
        % Go through all baseline days and compile stats
        FFvals_baseline=[];
        Tvals_baseline=[];
        for j=1:length(Params.SeqFilter.BaselineDays);
            day=Params.SeqFilter.BaselineDays(j);
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){day})
                continue
            end
            
            if ~isfield(AllDays_RegExpr.(['day_data' data_field]){day}, 'data_ParsedIntoSubclasses')
                continue;
            end
            
            ffvals=AllDays_RegExpr.(['day_data' data_field]){day}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals;
            tvals=AllDays_RegExpr.(['day_data' data_field]){day}.data_ParsedIntoSubclasses{i}.sub_class{ii}.Tvals;
            if isempty(ffvals);
            continue
            end
            
            ffvals=ffvals(:,1:length_subclass); % remove nans
            FFvals_baseline=[FFvals_baseline; ffvals];
            Tvals_baseline=[Tvals_baseline; tvals];
        end
        
        % =============== SLIDE IN BASELINE DATA INTO STRUCT
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals=FFvals_baseline;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.Tvals=Tvals_baseline;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.FFmean=mean(FFvals_baseline,1);
                
        % ============================ EACH DAY SUMMARY STATS - BOTH DAY
        % AND RELATIVE TO BASELINE
        for j=1:NumDays;
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j})
                continue
            end
            
            if ~isfield(AllDays_RegExpr.(['day_data' data_field]){j}, 'data_ParsedIntoSubclasses')
                continue;
            end
            
            
            % Today's data
            ffvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals;
            
            if isempty(ffvals);
                continue;
            end
            
            % remove nans
            ffvals=ffvals(:,1:length_subclass); 

            % baseline data
            ffval_baseline=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.FFmean;
            
            % ========================= SUMMARY STATS (not rel to baseline);
            numrows=size(ffvals,1);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFmean=mean(ffvals,1);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFsem=lt_sem(ffvals);
            
            % =========================== STATS Subtract Baseline
            if ~isempty(ffval_baseline);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base=ffvals-repmat(ffval_baseline,numrows,1);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_MEAN=mean(...
                AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_SEM=lt_sem(...
                AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base);
            else
                AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base=nan;
                AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_MEAN=nan;
                AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_SEM=nan;
            end
        end
        
    end
end


%% GET ONE SET OF STATS FOR EACH CLASS (COLLAPSING ALL SUBCLASSES INTO THE CLASS);
% Note: means will be aligned to first syllable (so only makes sense for
% open-ended repeats


for i=1:NumRegExprClasses;
    
    regexpr_class=Params.RegExpr.expressions{i};
    
    % ========= BASELINE Collect stats across all days
        FFvals_baseline=[];
        Tvals_baseline=[];
        for j=1:length(Params.SeqFilter.BaselineDays);
            day=Params.SeqFilter.BaselineDays(j);
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){day})
                continue
            end
            
             if ~isfield(AllDays_RegExpr.(['day_data' data_field]){day}.data_WithOutlier{i}, 'Final_ARRAYS')
                continue
            end
            
           
            ffvals=AllDays_RegExpr.(['day_data' data_field]){day}.data_WithOutlier{i}.Final_ARRAYS.FFvals;
            tvals=AllDays_RegExpr.(['day_data' data_field]){day}.data_WithOutlier{i}.Final_ARRAYS.Song_DateNum;
            
            if isempty(ffvals);
            continue
            end
            
            % === make the dim2 magnitude the same (i.e. pad the smaller matrix
            % with nan)
            if ~isempty(FFvals_baseline);
                if size(FFvals_baseline, 2)>size(ffvals,2);
                    % pad ffvals
                    tmp=size(FFvals_baseline, 2)-size(ffvals,2);
                    ffvals=[ffvals nan(size(ffvals,1),tmp)];
                elseif size(FFvals_baseline, 2)<size(ffvals,2);
                    % pad baseline
                    tmp=size(ffvals,2)-size(FFvals_baseline, 2);
                    FFvals_baseline=[FFvals_baseline nan(size(FFvals_baseline,1),tmp)];
                end
            end
                    
            % ==== Add to compiled data
            FFvals_baseline=[FFvals_baseline; ffvals];
            Tvals_baseline=[Tvals_baseline; tvals];
        end
        
        
        % =============== SLIDE IN BASELINE DATA INTO STRUCT
        AllDays_RegExpr.(['baseline' data_field]).data_WithOutlier{i}.FFvals=FFvals_baseline;
        AllDays_RegExpr.(['baseline' data_field]).data_WithOutlier{i}.Tvals=Tvals_baseline;
        AllDays_RegExpr.(['baseline' data_field]).data_WithOutlier{i}.FFmean=nanmean(FFvals_baseline,1);
        AllDays_RegExpr.(['baseline' data_field]).data_WithOutlier{i}.N=sum(~isnan(FFvals_baseline),1);
        
        % =================== EACH DAY SUMMARY STATS (day and relative to
        % baseline)
        
        for j=1:NumDays;
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j})
                continue
            end
            
            if ~isfield(AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}, 'Final_ARRAYS')
                continue
            end
            
            
            % Today's data
            ffvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.FFvals;
            
            if isempty(ffvals);
                continue;
            end
            
            % baseline data
            ffval_baseline=AllDays_RegExpr.(['baseline' data_field]).data_WithOutlier{i}.FFmean;
            
            % remove baseline columns that are not in today
            if size(ffval_baseline,2)>size(ffvals,2);
                tmp=size(ffval_baseline,2)-size(ffvals,2);
                ffval_baseline(end-tmp+1:end)=[];
            end
            
            % remove data columns that are not in baseline
            if size(ffval_baseline,2)<size(ffvals,2);
                tmp=size(ffvals,2)-size(ffval_baseline,2);
                ffvals(:,end-tmp+1:end)=[];
            end
            
            % ========================= SUMMARY STATS (not rel to baseline);
            numrows=size(ffvals,1);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.FFmean=nanmean(ffvals,1);

            N_mat=sum(~isnan(ffvals),1);
            tmp_std=nanstd(ffvals,0,1);
            tmp_sem=tmp_std./sqrt(N_mat-1);
            tmp_sem(N_mat<2)=nan;
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.FFsem=tmp_sem;
            
            
            % =========================== STATS Subtract Baseline
            if ~isempty(ffval_baseline);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base=ffvals-repmat(ffval_baseline,numrows,1);
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base_MEAN=...
                AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.FFmean-ffval_baseline;
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base_SEM=...
                AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.FFsem;
            else
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base=nan;
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base_MEAN=nan;
            AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base_SEM=nan;
            end
        end
        
end



%% PLOT ALL CLASSES (SUBCLASSES ALL COMBINED)
% RAW VALUES + BASELINE SUBTRACTED

% will collect day means in here:
FFmean_alldays_raw=cell(NumRegExprClasses, NumDays);
FFmean_alldays_minusbase=cell(NumRegExprClasses, NumDays);
FFsem_alldays_raw=cell(NumRegExprClasses, NumDays);
FFsem_alldays_minusbase=cell(NumRegExprClasses, NumDays);


for i=1:NumRegExprClasses;
    regexpr_string=Params.RegExpr.expressions{i};
    NumSylsMax=length(AllDays_RegExpr.(['baseline' data_field]).data_WithOutlier{i}.FFmean)+1;
    plot_cols=lt_make_plot_colors(NumSylsMax, 0, [1 0 0]);
    
    lt_figure; hold on;
    PlotLegend=0; % switches to 1 once legend is plotted;
        
    for j=1:NumDays;
        
        if isempty(AllDays_RegExpr.(['day_data' data_field]){j});
            continue;
        end
        
        
        % === GATHER DATA (RAW)
        FFvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.FFvals;
        Tvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.Song_DateNum;
        
        if isempty(FFvals) || size(FFvals,1)==1;
            continue;
        end
        
        % convert Tvals to days in expt
        eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
        Tvals=eventtimes.FinalValue;
        
        % get day means
        FFmean=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.FFmean;
        FFsem=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.FFsem;
        
        
        % GET RELATIVE TO BASELINE
        FFvals_minusB=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base;
        FFmean_minusB=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base_MEAN;
        FFsem_minusB=AllDays_RegExpr.(['day_data' data_field]){j}.data_WithOutlier{i}.Final_ARRAYS.STATS.REL_TO_BASELINE.FFvals_minus_base_SEM;
        
        
        % ==== COLLECT DAY MEANS AND SEMS TO LATER PLOT ACROSS DAYS (put
        % into a matrix, and fill in empty spots with nans
        FFmean_alldays_raw{i,j}=FFmean;
        FFmean_alldays_minusbase{i,j}=FFmean_minusB;
        
        FFsem_alldays_raw{i,j}=FFsem;
        FFsem_alldays_minusbase{i,j}=FFsem_minusB;
        
        % === PLOT
        % Parse into each syllable position, then plot
        numsyls=size(FFvals,2); % some motifs might be shorter
        hplot=[];
        for jj=1:numsyls;
            
            % ------- PLOT 1) PLOT ALL TRIALS
            lt_subplot(2,2,1); hold on;
            
            % remove nans
            inds=~isnan(FFvals(:, jj));
            
            X=Tvals(inds);
            Y=FFvals(inds, jj);
            if ~isempty(X);
                
                hplot(jj)=lt_plot(X, Y, {'Color',plot_cols{jj}});
                
                % Add day mean
                if length(X)>1;
                    try
                    errorbar(ceil(Tvals(1)), FFmean(jj), FFsem(jj), 's','MarkerSize',9,'Color',plot_cols{jj});
                    catch err
                    end
                end
                
                % ------- PLOT 2) PLOT REL TO BASELINE (vals)
                lt_subplot(2,2,2); hold on;
                try
                lt_plot(Tvals(inds), FFvals_minusB(inds,jj), {'Color',plot_cols{jj}});
                % add mean
                if length(X)>1;
                    errorbar(ceil(Tvals(1)), FFmean_minusB(jj), FFsem_minusB(jj), 's','MarkerSize',9,'Color',plot_cols{jj})
                end
                lt_plot_zeroline;
                catch err
                end
            end
        end
        
        if length(hplot)==NumSylsMax && PlotLegend==0; % otherwise legend will not have all syls
            try
                legend_cell={};
                for k=1:NumSylsMax;
                    legend_cell{k}=num2str(k);
                end
                legend(legend_cell);
                PlotLegend=1;
            catch err
            end
            
        end
        
    end
    
    % ===== PLOT MEANS ACROSS DAYS
    hplot=[];
    for jj=1:NumSylsMax;
        % === PLOT 3) RAW VALS
        % --------------------- collect data across days
        datacellarray_mean=FFmean_alldays_raw;
        datacellarray_sem=FFsem_alldays_raw;
        
        for jjj=1:NumDays;
                        
            TMP_mean=datacellarray_mean{i, jjj}; % vals for this day and reg exp string
            TMP_sem=datacellarray_sem{i, jjj};
            if length(TMP_mean)>=jj;
                tmp_mean=TMP_mean(jj);
                tmp_sem=TMP_sem(jj);
            else
                tmp_mean=nan;
                tmp_sem=nan;
            end
            Ymean(jjj)=tmp_mean;
            Ysem(jjj)=tmp_sem;
        end
        % ---------------------------------------------
        
        lt_subplot(2,2,3); hold on;
        shadedErrorBar(1:NumDays, Ymean, Ysem, {'o-','MarkerSize', 9,'MarkerFaceColor', plot_cols{jj}, 'Color',plot_cols{jj}},1);

        % === PLOT 4) MINUS BASELINE
        % --------------------- collect data across days
        datacellarray_mean=FFmean_alldays_minusbase;
        datacellarray_sem=FFsem_alldays_minusbase;
        
        for jjj=1:NumDays;
            
            TMP_mean=datacellarray_mean{i, jjj}; % vals for this day and reg exp string
            TMP_sem=datacellarray_sem{i, jjj};
            if length(TMP_mean)>=jj;
                tmp_mean=TMP_mean(jj);
                tmp_sem=TMP_sem(jj);                
            else
                tmp_mean=nan;
                tmp_sem=nan;                
            end
            Ymean(jjj)=tmp_mean;
            Ysem(jjj)=tmp_sem;
        end
        % ---------------------------------------------
        lt_subplot(2,2,4); hold on;
%         hplot(jj)=plot(1:NumDays, Ymean
        tmp=shadedErrorBar(1:NumDays, Ymean, Ysem, {'o-','MarkerSize', 9,'MarkerFaceColor', plot_cols{jj}, 'Color',plot_cols{jj}},1);
        hplot(jj)=tmp.mainLine;
        
        lt_plot_zeroline;

    end
    
    % legend
                    legend_cell={};
                for k=1:NumSylsMax;
                    legend_cell{k}=num2str(k);
                end

    legend(hplot, legend_cell)
    
    
    xlabel('days'); ylabel('FF (hz)');
    Fn_AnnotateWNLines(plotWNdays, ylim)
    
    lt_subtitle(regexpr_string);
end


%% PLOT SUBCLASSES PARSED - RAW DATA ACROSS DAYS - rewrite, it repeats 4 times.

% ==== PLOT ALL SUBCLASSES OF ALL CLASSES SEPARATELY
for i=1:NumRegExprClasses;
    regexpr_string=Params.RegExpr.expressions{i};
    
    num_subclasses=length(Params.RegExpr.subexpressions{i});

    if num_subclasses<2; % then skip this analysis - is redundant with class analysis
        continue;
    end
    
    for ii=1:num_subclasses;
        subclass=Params.RegExpr.subexpressions{i}{ii};
        
        % === 1) PLOT learning across days (actual FF)
        lt_figure; hold on;
        
        % == Raw values + means;
        lt_subplot(2,2,1); hold on;
        title(subclass);
        
        for j=1:NumDays;
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j});
                continue;
            end
            
            
            Tvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.Tvals;
            FFvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals;
            
            if isempty(FFvals) || size(FFvals,1)==1;
                continue;
            end
            
            % convert Tvals to days in expt
            eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
            Tvals=eventtimes.FinalValue;
            
            % get day means
            FFmean=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFmean;
            FFsem=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFsem;
            
            % Parse into each syllable, then plot
            numsyls=length(subclass);
            plot_cols=lt_make_plot_colors(numsyls, 0, [1 0 0]);
            hplot=[];
            for jj=1:numsyls;
                
                % ======== PLOT ALL TRIALS
                lt_plot(Tvals, FFvals(:,jj), {'Color',plot_cols{jj}});
                
                % ======== PLOT DAY MEAN
                hplot(jj)=errorbar(ceil(Tvals(1)), FFmean(jj), FFsem(jj), 's','MarkerSize',9,'Color',plot_cols{jj});
                
            end
            
        end
        
        try
        legend_cell={};
        for k=1:length(subclass);
            legend_cell{k}=subclass(k);
        end
        legend(hplot, legend_cell);
        catch err
        end

        xlabel('days'); ylabel('FF (hz)');
        Fn_AnnotateWNLines(plotWNdays, ylim)
        
        % ==== MEAN ONLY
        lt_subplot(2,2,2); hold on;
        FFmean_matrix=[];
        Tval_matrix=[];
        FFsem_matrix=[];
        for j=1:NumDays;
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j});
                continue;
            end
            
            Tvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.Tvals;
            FFvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals;
            
            if isempty(FFvals) || size(FFvals,1)==1;
                continue;
            end
            
            % convert Tvals to days in expt
            eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
            Tvals=eventtimes.FinalValue;
            
            % get day means
            FFmean=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFmean;
            FFsem=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFsem;
            
            % ==== Across days matrix
            FFmean_matrix=[FFmean_matrix; FFmean];
            FFsem_matrix=[FFsem_matrix; FFsem];
            Tval_matrix=[Tval_matrix; j];
        end
        
        if isempty(FFmean_matrix);
            continue
        end
        
        % Parse into each syllable, then plot
        numsyls=length(subclass);
        plot_cols=lt_make_plot_colors(numsyls, 0, [1 0 0]);
        
        for jj=1:numsyls;
            
            % ======== PLOT DAY MEAN
            if length(Tval_matrix)>1
            shadedErrorBar(Tval_matrix, FFmean_matrix(:,jj), FFsem_matrix(:, jj), {'o','MarkerFaceColor',plot_cols{jj},'MarkerSize',9,'Color',plot_cols{jj}},1);
            end
        end
        
        xlabel('days'); ylabel('FF (hz)');
        Fn_AnnotateWNLines(plotWNdays, ylim)

        
        
        % ============================================= 2) PLOT REL BASELINE
        % == Raw values + means;
        lt_subplot(2,2,3); hold on;

        for j=1:NumDays;
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j});
                continue;
            end
            if ~isfield(AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii},'REL_TO_BASELINE');
                continue;
            end
           
            
            Tvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.Tvals;
            FFvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base;
            
            if isempty(FFvals) || size(FFvals,1)==1;
                continue;
            end
            
            % convert Tvals to days in expt
            eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
            Tvals=eventtimes.FinalValue;
            
            % get day means
            FFmean=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_MEAN;
            FFsem=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_SEM;
            
            % Parse into each syllable, then plot
            numsyls=length(subclass);
            plot_cols=lt_make_plot_colors(numsyls, 0, [1 0 0]);
            hplot=[];
            for jj=1:numsyls;
                
                % ======== PLOT ALL TRIALS
                lt_plot(Tvals, FFvals(:,jj), {'Color',plot_cols{jj}});
                
                % ======== PLOT DAY MEAN
            hplot(jj)=errorbar(j+0.5, FFmean(jj), FFsem(jj), 's','MarkerSize',9,'Color',plot_cols{jj});
                
            end
            
        end
        
        try
        legend_cell={};
        for k=1:length(subclass);
            legend_cell{k}=subclass(k);
        end
        legend(hplot, legend_cell);
        catch err
        end
        xlabel('days'); ylabel('FF (hz)');
        Fn_AnnotateWNLines(plotWNdays, ylim)
line(xlim, [0 0]);

        
        % ==== MEAN ONLY
        lt_subplot(2,2,4); hold on;
        FFmean_matrix=[];
        Tval_matrix=[];
        FFsem_matrix=[];
        for j=1:NumDays;
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j});
                continue;
            end
            if ~isfield(AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii},'REL_TO_BASELINE');
                continue;
            end

            Tvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.Tvals;
            FFvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base;
            
            if isempty(FFvals) || size(FFvals,1)==1;
                continue;
            end
            
            % convert Tvals to days in expt
            eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
            Tvals=eventtimes.FinalValue;
            
            % get day means
            FFmean=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_MEAN;
            FFsem=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_SEM;
            
            % ==== Across days matrix
            FFmean_matrix=[FFmean_matrix; FFmean];
            FFsem_matrix=[FFsem_matrix; FFsem];
            Tval_matrix=[Tval_matrix; j];
        end
        
        if isempty(FFmean_matrix);
            continue
        end
        
        % Parse into each syllable, then plot
        numsyls=length(subclass);
        plot_cols=lt_make_plot_colors(numsyls, 0, [1 0 0]);
        
        for jj=1:numsyls;
            % ======== PLOT DAY MEAN
            if length(Tval_matrix)>1
            shadedErrorBar(Tval_matrix, FFmean_matrix(:,jj), FFsem_matrix(:, jj), {'o','MarkerFaceColor',plot_cols{jj},'MarkerSize',9,'Color',plot_cols{jj}},1);
            end
        end
                xlabel('days'); ylabel('FF (hz)');
        Fn_AnnotateWNLines(plotWNdays, ylim)
line(xlim, [0 0]);

    end
end

            
%% PLOT ON ONE FIGURE ALL SUBCLASSES

for i=1:NumRegExprClasses;
    regexpr_string=Params.RegExpr.expressions{i};
    num_subclasses=length(Params.RegExpr.subexpressions{i});
    
    if num_subclasses<2;
        continue; 
    end
        
    lt_figure; hold on;
    
    hsb=[];
    for ii=1:num_subclasses;
        
        subclass=Params.RegExpr.subexpressions{i}{ii};
        plot_cols=lt_make_plot_colors(length(Params.RegExpr.subexpressions{i}{end}), 0, [1 0 0]);
        
        % === PLOT DAY MEANS (MINUS BASELINE)
        FFmean_matrix=[];
        Tval_matrix=[];
        FFsem_matrix=[];
        for j=1:NumDays;
            
            if isempty(AllDays_RegExpr.(['day_data' data_field]){j});
                continue;
            end
            if ~isfield(AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii},'REL_TO_BASELINE');
                continue;
            end
            
            FFvals=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base;
            
            if isempty(FFvals) || size(FFvals,1)==1;
                continue;
            end
            
            % get day means
            FFmean=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_MEAN;
            FFsem=AllDays_RegExpr.(['day_data' data_field]){j}.data_ParsedIntoSubclasses{i}.sub_class{ii}.REL_TO_BASELINE.FFvals_minus_base_SEM;
            
            % ==== Across days matrix
            FFmean_matrix=[FFmean_matrix; FFmean];
            FFsem_matrix=[FFsem_matrix; FFsem];
            Tval_matrix=[Tval_matrix; j];
        end
        
        if isempty(FFmean_matrix);
            continue
        end
        
        % Parse into each syllable, then plot
        numsyls=length(subclass);
        
        hplot=[];
        hsb(ii)=lt_subplot(ceil(num_subclasses/2), 2, ii); hold on;
        title(subclass);
        
        for jj=1:numsyls;
            hplot(jj)= plot(Tval_matrix, FFmean_matrix(:,jj), 'o', 'Color', plot_cols{jj}, 'MarkerFaceColor',plot_cols{jj}); % only to get legend
            if length(Tval_matrix)>1;
                % ======== PLOT DAY MEAN
                shadedErrorBar(Tval_matrix, FFmean_matrix(:,jj), FFsem_matrix(:, jj), {'o','MarkerFaceColor',plot_cols{jj},'MarkerSize',9,'Color',plot_cols{jj}},1);
            end
        end
        
        % legend
        if ii==num_subclasses;
                legend_cell={};
                for k=1:length(subclass);
                    legend_cell{k}=subclass(k);
                end
                legend(hplot, legend_cell);
        end
        
        xlabel('days'); ylabel('FF (hz)');
        Fn_AnnotateWNLines(plotWNdays, ylim)
        line(xlim, [0 0], 'Color','k');
        
    end
    
    linkaxes(hsb, 'xy');
end


%% ONE FIGURE FOR ALL SUBCLASSES, BUT BARPLOT FOR END OF LEARNING VALUES
% IN PROGRESS
           

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++ CORRELATIONS

%% CORRELATIONS, using mean pitch in the song (to be able to compare syls in different motifs)

day_data_field=['day_data' data_field];
baseline_field=['baseline' data_field];

% ================== Extracting data - each row is one song, and each column is one type of
% syl.

% method: go through each song one by one. for each song find all classes
% and subclasses that have data for this song. collect their data.

% 1) Collect all songs
AllSongs_baseline=[];
baseline_days=Params.SeqFilter.BaselineDays;
for i=1:length(baseline_days);
    day=baseline_days(i); 
    
    if isempty(AllDays_RegExpr.(['day_data' data_field]){day});
       continue;
    end
    
    if size(AllDays_RegExpr.(['day_data' data_field]){day}.AllSong_datenums,1)==1; % then it is in columns, want it in one column

        AllSongs_baseline=[AllSongs_baseline; AllDays_RegExpr.(['day_data' data_field]){day}.AllSong_datenums'];
    else
        AllSongs_baseline=[AllSongs_baseline; AllDays_RegExpr.(['day_data' data_field]){day}.AllSong_datenums];
    end
        
end

AllSongs_baseline=unique(AllSongs_baseline);


% 2) For each song, go through all classes and subclasses and collect data
% -- First, initiate the output matrices
for j=1:NumRegExprClasses;
    % ==== CLASSES (subclasses compressed)
    tmp=AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.FFvals;
    AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd=nan(length(AllSongs_baseline), size(tmp,2)); % prepare the output ffmean val matrix
    
    % ==== SUBCLASSES
    num_subclasses = length(AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class);
    for jj=1:num_subclasses;
        tmp=AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.FFvals;
        AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd=nan(length(AllSongs_baseline), size(tmp,2)); % prepare the output ffmean val matrix
    end
end


for i=1:length(AllSongs_baseline);
    song_datenum=AllSongs_baseline(i);
    
    % ================ GO THROUGH ALL CLASSES (compressed, i.e. no
    % subclasses)
    for j=1:NumRegExprClasses;
        
        curr_song_Inds= find(AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.Tvals==song_datenum);
        ffvals=AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.FFvals(curr_song_Inds, :); % for all syls in this class
        
        
        % --- OUTPUT - All Raw Data
        AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.CORRELATIONS.SONG_BY_SONG.FFvals_raw{i}=ffvals;
        AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.CORRELATIONS.SONG_BY_SONG.IndsOfSong_raw{i}=curr_song_Inds;
        AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.CORRELATIONS.SONG_BY_SONG.Song_datenum(i)=song_datenum;
        
        % --- OUTPUT - take mean for this song and put in matrix
        if size(ffvals,1)>1; % if multiple rends in this song
            ffmean=nanmean(ffvals);
        elseif size(ffvals, 1)==0; % then no rendition of this class in this song - put nans
            ffmean=nan(1, size(ffvals,2));
        else ffmean=ffvals; % then has exactly one rendition, don't take mean
        end
        
        AllDays_RegExpr.(baseline_field).data_WithOutlier{j}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd(i, :) = ffmean;
    end
    
    % ================ GO THROUGH ALL SUBCLASSES (for each CLASS
    for j=1:NumRegExprClasses;
        num_subclasses = length(AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class);
        
        for jj=1:num_subclasses;
          
            curr_song_Inds=find(AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.Tvals==song_datenum);
            ffvals=AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.FFvals(curr_song_Inds, :);
            
            
            % --- OUTPUT - All Raw Data
            AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.CORRELATIONS.SONG_BY_SONG.FFvals_raw{i}=ffvals;
            AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.CORRELATIONS.SONG_BY_SONG.IndsOfSong_raw{i}=curr_song_Inds;
            AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.CORRELATIONS.SONG_BY_SONG.Song_datenum{i}=song_datenum;
            
            % --- OUTPUT - Take mean for this song and put into matrix
            % (global song ind)
            if size(ffvals,1)>1; % if multiple rends in this song
                ffmean=nanmean(ffvals);
            elseif size(ffvals, 1)==0; % then no rendition of this class in this song - put nans
                ffmean=nan(1, size(ffvals,2));
            else ffmean=ffvals; % then has exactly one rendition, don't take mean
            end
            
            AllDays_RegExpr.(baseline_field).data_ParsedIntoSubclasses{j}.sub_class{jj}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd(i, :) = ffmean;
            
        end
    end
end

        
%% CALCUALATE AND PLOT CORRELATIONS USING SONG MEANS (majore classes, ignoreing subclasses)

% ===== COMPUTE GLOBAL FF matrix (each row is one global song ind, so some
% nan are possible)
SylsInOrder_Compiled={};
FFmean_by_song_ind_Compiled=[];

for i=1:NumRegExprClasses;
    
    % Collect data
    SylsInOrder=Params.RegExpr.expressions{i};
    FFmean_by_song_ind=AllDays_RegExpr.(baseline_field).data_WithOutlier{i}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd;
    
    % Put into output matrix and cell
    for ii=1:length(SylsInOrder);
        SylsInOrder_Compiled=[SylsInOrder_Compiled SylsInOrder(ii)];
    end
    FFmean_by_song_ind_Compiled=[FFmean_by_song_ind_Compiled FFmean_by_song_ind];
end

%  ==== REMOVE ANY SONGS THAT HAVE NAN FOR ANY OF THE SYLS
% IGNORE ANY COLUMNS THAT ARE ALL nans - those are just syls that don't
% have ffval data.

% what columns have all nan?
cols_to_ignore=find(sum(isnan(FFmean_by_song_ind_Compiled),1)==size(FFmean_by_song_ind_Compiled,1));

% remove rows that have any nan(other than column to ignore)
TestMatrix=FFmean_by_song_ind_Compiled;
TestMatrix(:,cols_to_ignore)=[];

Inds_nan=isnan(TestMatrix);
Inds_nan=sum(Inds_nan,2)>0; % get the rows with any nan

% remove those inds
disp(['Removed ' num2str(sum(Inds_nan)) ' rows (songs) out of ' num2str(length(Inds_nan)) ' total, as missed at least one syl']);

FFmean_by_song_ind_Compiled_NanSongsRemoved=FFmean_by_song_ind_Compiled;
FFmean_by_song_ind_Compiled_NanSongsRemoved(Inds_nan,:)=[];

% ===== CALCULATE AND PLOT
lt_figure; hold on;
title('Correlation matrix between all classes (using mean in song)');
[RhoMat, PvalMat]=corr(FFmean_by_song_ind_Compiled_NanSongsRemoved);

% --- OUTPUT STRUCTURE
        AllDays_RegExpr.(baseline_field).Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled=FFmean_by_song_ind_Compiled;
        AllDays_RegExpr.(baseline_field).Correlations_Across_Classes.SONG_BY_SONG.Song_datenums=...
            AllDays_RegExpr.(baseline_field).data_WithOutlier{i}.CORRELATIONS.SONG_BY_SONG.Song_datenum;
        
        AllDays_RegExpr.(baseline_field).Correlations_Across_Classes.SONG_BY_SONG.FFmean_by_song_ind_Compiled_NanSongsRemoved=FFmean_by_song_ind_Compiled_NanSongsRemoved;
        
        AllDays_RegExpr.(baseline_field).Correlations_Across_Classes.SONG_BY_SONG.SylsInOrder_Compiled=SylsInOrder_Compiled;
        
        AllDays_RegExpr.(baseline_field).Correlations_Across_Classes.SONG_BY_SONG.RhoMat=RhoMat;
        AllDays_RegExpr.(baseline_field).Correlations_Across_Classes.SONG_BY_SONG.PvalMat=PvalMat;
        
        
% ==========================================================
% === PLOT
% 1) Rho
colormap('gray');
imagesc(RhoMat,[-0.15 0.7]);
set(gca,'YTick',1:length(SylsInOrder_Compiled), 'YTickLabel', SylsInOrder_Compiled);
set(gca,'XTick',1:length(SylsInOrder_Compiled), 'XTickLabel', SylsInOrder_Compiled);

colorbar;


% ============================== PLOT AS BAR PLOTS
lt_figure; hold on;

for i=1:length(SylsInOrder_Compiled);
    lt_subplot(ceil(length(SylsInOrder_Compiled)/4), 4, i); hold on;
    title(['syl: ' SylsInOrder_Compiled{i}])
    
    Y=RhoMat(:,i);
    
    bar(Y);
    xlim([0 length(SylsInOrder_Compiled)+1]);
    ylim([-0.15 0.7])
    set(gca,'XTick',1:length(SylsInOrder_Compiled), 'XTickLabel', SylsInOrder_Compiled);
end
lt_subtitle('Correlation between all classes (song mean)');


        
%% CORRELATIONS USING MOTIFS (baseline), all subclasses? With and Without subtracting song mean
MinRendPerSong=3; % min rends of motif to take data for song-subtracted correaltions.

% ==== COLLECT CORRELATION MATRICES
for i=1:NumRegExprClasses;
            
    regexpr_string=Params.RegExpr.expressions{i};
    
    num_subclasses=length(Params.RegExpr.subexpressions{i});
    
    % For each subclass:
    for ii=1:num_subclasses;
        
        subclass=Params.RegExpr.subexpressions{i}{ii};
        
        % ============================================= EXTRACT DATA ===========================================================
        % ++++++++++++++++++++++ DATA MATRIX:
        FFmatrix_base=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals;
        
        % === Make sure enough data exists (at least 8 trials)
        if size(FFmatrix_base,1)<9;
            continue;
        end
        
        
        % +++++++++++++++++++++++ CORRELATIONS (ACROSS MOTIFS)
        % ========= GET CORRELATION MATRIX
        [RhoMat, PvalMat]=corr(FFmatrix_base);
        
        % --- OUTPUT STRUCTURE
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.RhoMat=RhoMat;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.PvalMat=PvalMat;
        
        
        % +++++++++++++++++++++++++ CORRELATIONS (SUBTRACTING SONG MEAN
        % FIRST)
        % 1) collect all raw values, but before compiling subtract mean of
        % the song
        num_songs=size(AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd,1);
        ffvals_MinMean_compiled=[];
        for j=1:num_songs;
            ffvals=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SONG_BY_SONG.FFvals_raw{j};
            
            if isempty(ffvals) | any(isnan(ffvals));
                % then no data in this song
                continue
            end
            
            ffvals_mean=mean(ffvals,1);
            ffvals_minusmean=ffvals-repmat(ffvals_mean, size(ffvals,1), 1);
            
            ffvals_MinMean_compiled=[ffvals_MinMean_compiled; ffvals_minusmean];
        end
        
        % --- CORRELATIONS
        [RhoMat_MinSongMean, PvalMat_MinSongMean]=corr(ffvals_MinMean_compiled);

        
        % --- OUTPUT STRUCTURE
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN.ffvals=ffvals_MinMean_compiled;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN.RhoMat=RhoMat_MinSongMean;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN.PvalMat=PvalMat_MinSongMean;
        
        
        % +++++++++++++++++++++++++ CORRELATIONS (SUBTRACTING SONG MEAN
        % FIRST, ONLY TAKING IF SONG HAS OVER N RENDITIONS OF MOTIF)
        % 1) collect all raw values, but before compiling subtract mean of
        % the song
        num_songs=size(AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SONG_BY_SONG.FFmean_matrix_by_GlobalSongInd,1);
        ffvals_MinMean_NoShortSong_compiled=[];
        songs_taken_compiled=[];
        for j=1:num_songs;
            ffvals=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SONG_BY_SONG.FFvals_raw{j};
            
            if isempty(ffvals) | any(isnan(ffvals));
                % then no data in this song
                continue
            end
            
            if size(ffvals,1)<MinRendPerSong;
                % then skip this song, as does not have enough rends of
                % motif
                continue;
            end
             
            ffvals_mean=mean(ffvals,1);
            ffvals_minusmean=ffvals-repmat(ffvals_mean, size(ffvals,1), 1);
            
            ffvals_MinMean_NoShortSong_compiled=[ffvals_MinMean_NoShortSong_compiled; ffvals_minusmean];
            songs_taken_compiled=[songs_taken_compiled j];
        end
        
        % --- CORRELATIONS
        if ~isempty(ffvals_MinMean_NoShortSong_compiled);
        [RhoMat_MinSongMean_NoShortSong, PvalMat_MinSongMean_NoShortSong]=corr(ffvals_MinMean_NoShortSong_compiled);
        end
        
        % --- OUTPUT STRUCTURE
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.ffvals=ffvals_MinMean_NoShortSong_compiled;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.Songs_kept=songs_taken_compiled;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.RhoMat=RhoMat_MinSongMean_NoShortSong;
        AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.PvalMat=PvalMat_MinSongMean_NoShortSong;

        % ========================================== PLOT =====================================================================
        % +++++++++++++++++++++++ 1) WITHOUT SUBTRACTING SONG MEANS 
        lt_figure; hold on;

        % 1) Rho
        lt_subplot(2, num_subclasses, ii); hold on;
        title(['Rho: ' subclass ' (N= ' num2str(size(FFmatrix_base,1)) ')']);
        
        colormap('gray');
        imagesc(RhoMat,[-0.15 0.7]);
        set(gca,'YTick',1:length(subclass));
        %         set(gca,'YTickLabel',subclass);
        set(gca,'XTick',1:length(subclass));
        %         set(gca,'XTickLabel',subclass);
        %         b = strread(subclass,'%s')
        % b = sscanf(subclass,'%c')
        
        colorbar;
        
        % 2) P-value
        lt_subplot(2, num_subclasses, ii+num_subclasses); hold on;
        title(['P-value: ' subclass ' (N= ' num2str(size(FFmatrix_base,1)) ')']);
        
        colormap('gray');
        imagesc(PvalMat,[0, 0.1]);
        set(gca,'YTick',1:length(subclass));
        %         set(gca,'YTickLabel',{subclass});
        set(gca,'XTick',1:length(subclass));
        %         set(gca,'XTickLabel',subclass);
        colorbar;
        
        lt_subtitle(['Correlations (using motifs): ' regexpr_string]);

        
        % +++++++++++++++++++++++ 2) AFTER SUBTRACTING SONG MEANS 
        lt_figure; hold on;

        % 1) Rho
        lt_subplot(2, num_subclasses, ii); hold on;
        title(['Rho: ' subclass ' (N= ' num2str(size(FFmatrix_base,1)) ')']);
        
        colormap('gray');
        imagesc(RhoMat_MinSongMean,[-0.15 0.7]);
        set(gca,'YTick',1:length(subclass));
        %         set(gca,'YTickLabel',subclass);
        set(gca,'XTick',1:length(subclass));
        %         set(gca,'XTickLabel',subclass);
        %         b = strread(subclass,'%s')
        % b = sscanf(subclass,'%c')
        
        colorbar;
        
        % 2) P-value
        lt_subplot(2, num_subclasses, ii+num_subclasses); hold on;
        title(['P-value: ' subclass ' (N= ' num2str(size(FFmatrix_base,1)) ')']);
        
        colormap('gray');
        imagesc(PvalMat_MinSongMean,[0, 0.1]);
        set(gca,'YTick',1:length(subclass));
        %         set(gca,'YTickLabel',{subclass});
        set(gca,'XTick',1:length(subclass));
        %         set(gca,'XTickLabel',subclass);
        colorbar;
        
    lt_subtitle(['Correlations (using motifs, song mean subtracted): ' regexpr_string]);

            % +++++++++++++++++++++++ 3) AFTER SUBTRACTING SONG MEANS (No
            % Short Songs)
        lt_figure; hold on;

        % 1) Rho
        lt_subplot(2, num_subclasses, ii); hold on;
        title(['Rho: ' subclass ' (N= ' num2str(size(FFmatrix_base,1)) ')']);
        
        colormap('gray');
        imagesc(RhoMat_MinSongMean_NoShortSong,[-0.15 0.7]);
        set(gca,'YTick',1:length(subclass));
        %         set(gca,'YTickLabel',subclass);
        set(gca,'XTick',1:length(subclass));
        %         set(gca,'XTickLabel',subclass);
        %         b = strread(subclass,'%s')
        % b = sscanf(subclass,'%c')
        
        colorbar;
        
        % 2) P-value
        lt_subplot(2, num_subclasses, ii+num_subclasses); hold on;
        title(['P-value: ' subclass ' (N= ' num2str(size(FFmatrix_base,1)) ')']);
        
        colormap('gray');
        imagesc(PvalMat_MinSongMean_NoShortSong,[0, 0.1]);
        set(gca,'YTick',1:length(subclass));
        %         set(gca,'YTickLabel',{subclass});
        set(gca,'XTick',1:length(subclass));
        %         set(gca,'XTickLabel',subclass);
        colorbar;
        
    lt_subtitle(['Correlations (motifs, minus song mean, no short song): ' regexpr_string]);
    end
   
end

%% ========================================================
%% PICK A CLASS AND SUBCLASS AND PLOT THINGS IN MORE DETAILS
disp(' ');
disp('Will perform further analyses, plots, on selected regexp class and subclass');
disp(' ');

% ==== PERFORM ANALYSIS ONCE FOR EACH CLASS. If there are multiple
% subclasses in a class, then have the user choose.

for n=1:NumRegExprClasses;
    subclass_choice=[];
    
    % ========================= CLASS CHOICE:
    subclass_choice(1)=n;
    
    % ====================== SUBCLASS CHOICE
    num_subclasses=length(Params.RegExpr.subexpressions{subclass_choice(1)});
    
    if num_subclasses==1;
        % ==== AUTO CHOOSE (only one subclass)
        subclass_choice(2)=1;
    else
        % multiple subclasses present
        if SuppressQueries==0;
            % ======= USER CHOOSE (if there are multiple subclasses)
            % disp(Params.RegExpr.expressions)
            %
            % subclass_choice(1)=input('which class do you choose? (e.g. 1, 2, ...) ');
            %
            disp(Params.RegExpr.subexpressions{subclass_choice(1)})
            subclass_choice(2)=input('which subclass do you choose? (e.g. 1, 2, ...) ');
        else
            % ===== SUPPRESS USER CHOICE, Choose the middle subclass
            subclass_choice(2)=ceil(length(Params.RegExpr.subexpressions{subclass_choice(1)})/2);
        end
    end
        
        disp(' ');
        disp(['Current class: ' Params.RegExpr.expressions{subclass_choice(1)}]);
    disp(['Current subclass: ' Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)}]);
    pause(1);
    
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % ++++++++++++++++++++++++++++++++++++++++++++++++++ PLOT BAR PLOTS
    % 1) WITHOUT SUBTRACTING SONG MEAN
    num_syls=length(Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)});
    syllist=Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)};
    
    lt_figure; hold on;
    for i=1:num_syls;
        lt_subplot(ceil(num_syls/3), 3, i); hold on;
        title(['syl: ' syllist(i)]);
        
        Y=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.CORRELATIONS.RhoMat(i,:);
        
        bar(Y);
        ylim([-0.15 0.7]);
        xlim([0.5 num_syls+0.5]);
        set(gca, 'XTick', 1:num_syls);
    end
    
    lt_subtitle([Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)} '; rho (motifs)']);
    
    % 2) WITH SUBTRACTING SONG MEAN
    lt_figure; hold on;
    for i=1:num_syls;
        lt_subplot(ceil(num_syls/3), 3, i); hold on;
        title(['syl: ' syllist(i)]);
        
        Y=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.CORRELATIONS.SUBTRACT_SONG_MEAN.RhoMat(i,:);
        
        bar(Y);
        ylim([-0.15 0.7]);
        xlim([0.5 num_syls+0.5]);
        set(gca, 'XTick', 1:num_syls);
    end
    
    lt_subtitle([Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)} '; rho (song mean subtracted)']);
    
    % 2) WITH SUBTRACTING SONG MEAN
    lt_figure; hold on;
    for i=1:num_syls;
        lt_subplot(ceil(num_syls/3), 3, i); hold on;
        title(['syl: ' syllist(i)]);
        
        Y=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.RhoMat(i,:);
        
        bar(Y);
        ylim([-0.15 0.7]);
        xlim([0.5 num_syls+0.5]);
        set(gca, 'XTick', 1:num_syls);
    end
    
    lt_subtitle([Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)} '; rho (minus song mean, no short songs)']);
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % +++++++++++++++++++++++++++++++++++++++++++= PLOT ACROSS TIME
    
    if (0); % takes a long time, ignore for now.
    % === EXTRACT DATASET AND STATS
    regexpr_string=Params.RegExpr.expressions{subclass_choice(1)};
    subclass=Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)};
    
    % DATA MATRIX:
    FFmatrix_base=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.FFvals;
    Tmatrix_base=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.Tvals;
    num_syls=size(FFmatrix_base,2);
    plotcols=lt_make_plot_colors(num_syls,0,0);
    
    % === PLOT - each rendition the same color for all syls.
    lt_figure; hold on;
    counter=1;
    hsplot=[];
    
    numrends=length(Tmatrix_base);
    
    for i=1:numrends;
        
        color=[rand rand rand];
        
        for kk=1:num_syls;
            hsplot(kk)=lt_subplot(ceil(num_syls/3),3,kk); hold on;
            
            ffval=FFmatrix_base(i, kk);
            tval=Tmatrix_base(i);
            
            plot(tval, ffval, 'o', 'Color', color, 'MarkerFaceColor', color);
            
        end
    end
    
    % plot line for means
    for kk=1:num_syls;
        hsplot(kk)=lt_subplot(ceil(num_syls/3),3,kk); hold on;
        
        ffvals=FFmatrix_base(:,kk);
        
        % line for average ff
        ffmean=mean(ffvals);
        line(xlim, [ffmean ffmean], 'Color' ,'k');
        
    end
    
    lt_subtitle([subclass '; all rends diff color, all same motif rends same color']);
    linkaxes(hsplot, 'x');
    end
    
    % ++++++++++++++++++++++++++++++++++++++++++++++++
    % ++++++++++++++++++++++++= PLOT SCATTER - one syl against other syls
    % TO DO - subtract out drift, etc.  do I really see correlation between all
    % syls?
    
    % === EXTRACT DATASET AND STATS
    regexpr_string=Params.RegExpr.expressions{subclass_choice(1)};
    subclass=Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)};
    
    % DATA MATRIX:
    FFmatrix_base=AllDays_RegExpr.(['baseline' data_field]).data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.FFvals;
    num_syls=size(FFmatrix_base,2);
    plotcols=lt_make_plot_colors(num_syls,0,0);
    
    
    % PLOT each syl as a target syl, and scatter all the others
    % against it.
    lt_figure; hold on;
    counter=1;
    for kk=1:num_syls;
        
        for kkk=1:num_syls;
            
            lt_subplot(num_syls,num_syls,counter); hold on;
            
            targsyl_ff=FFmatrix_base(:,kk);
            othersyl_ff=FFmatrix_base(:,kkk);
            
            [b,bint,r,rint,stats,SummaryStats]=lt_regress(othersyl_ff,targsyl_ff,1,0);
            
            %         plot(targsyl_ff, othersyl_ff, 'ok', 'MarkerSize',4);
            %         plot(xlim,b(1) + b(2).*xlim,'-r','LineWidth',2);
            
            title(['syl ' num2str(kkk) ' vs. syl ' num2str(kk)]);
            %         ylabel(    num2str(kkk));
            %         xlabel(num2str(kk))
            counter=counter+1;
        end
    end
    
    lt_subtitle(subclass)
    
end

% OLD VERSION - went thru all subclasses and classes - exploratory
% for i=1:NumRegExprClasses;
%     
%     %     lt_figure; hold on;
%     regexpr_string=Params.RegExpr.expressions{i};
%     
%     num_subclasses=length(Params.RegExpr.subexpressions{i});
%     
%     % For each subclass:
%     for ii=1:num_subclasses;
%         
%         subclass=Params.RegExpr.subexpressions{i}{ii};
%         % DATA MATRIX:
%         FFmatrix_base=AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{i}.sub_class{ii}.FFvals;
%         
%         % === Make sure enough data exists (at least 8 trials)
%         if size(FFmatrix_base,1)<9;
%             continue;
%         end
%         
%         % ========= GET CORRELATION MATRIX
%         [RhoMat, PvalMat]=corr(FFmatrix_base);
%         
%         % --- OUTPUT STRUCTURE
%         AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.CorrMat=RhoMat;
%         AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{i}.sub_class{ii}.CORRELATIONS.CorrMat=PvalMat;
%         
%         
%         % ==========================================================
%         % === PLOT
%         % 0) Scatter, each syl as a target syl, and scatter all the others
%         % against it.
%         num_syls=size(FFmatrix_base,2);
%         plotcols=lt_make_plot_colors(num_syls,0,0);
%         
%         lt_figure; hold on;
%         for kk=1:num_syls;
%             lt_subplot(ceil(num_syls/3),3,kk); hold on;
%             title(['all against syl ' num2str(kk)]);
%             
%             for kkk=1:num_syls;
%                 if kkk==kk;
%                     continue
%                 end
%                 
%                 targsyl_ff=FFmatrix_base(:,kk);
%                 othersyl_ff=FFmatrix_base(:,kkk);
%                 
%                 lt_plot(targsyl_ff, othersyl_ff, {'Color',plotcols{kkk}});
%             end
%         end
%         
%         lt_subtitle(subclass)
%         
%         
%     end
% end



  
%% SCATTER PLOT, but for each rendition first subtract the mean of that syl in the song
% IN PROGRESS - first o code to get mean for each song. plot that and get
% corr for that. then do song-mean-subtracted.

% % === EXTRACT DATASET AND STATS
% regexpr_string=Params.RegExpr.expressions{subclass_choice(1)};
% subclass=Params.RegExpr.subexpressions{subclass_choice(1)}{subclass_choice(2)};
% 
% % DATA MATRIX:
% FFmatrix_base=AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{subclass_choice(1)}.sub_class{subclass_choice(2)}.FFvals;
% num_syls=size(FFmatrix_base,2);
% plotcols=lt_make_plot_colors(num_syls,0,0);
% 
% % ==========================================================
% % === PLOT
% % 0) Scatter, each syl as a target syl, and scatter all the others
% % against it.
% 
% lt_figure; hold on;
% counter=1;
% for kk=1:num_syls;
%     
%     for kkk=1:num_syls;
%         
%         lt_subplot(num_syls,num_syls,counter); hold on;
%         
%         targsyl_ff=FFmatrix_base(:,kk);
%         othersyl_ff=FFmatrix_base(:,kkk);
%         
%         [b,bint,r,rint,stats,SummaryStats]=lt_regress(othersyl_ff,targsyl_ff,1,0);
%         
% %         plot(targsyl_ff, othersyl_ff, 'ok', 'MarkerSize',4);
% %         plot(xlim,b(1) + b(2).*xlim,'-r','LineWidth',2);
% 
%         title(['syl ' num2str(kkk) ' vs. syl ' num2str(kk)]);
% %         ylabel(    num2str(kkk));
% %         xlabel(num2str(kk))
%         counter=counter+1;
%     end
% end
% 
% lt_subtitle(subclass)

%% ======= CORRELATE PITCH OF SYL IN ONE MOTIF RENDITION WITH OTHER MOTIF RENDITIONS IN THE SAME SONG
% Perform once for each syllable (defined by position in sequence)
% IN PROGRESS


% for i=1:NumRegExprClasses;
%     
%     numsyls_max=length(AllDays_RegExpr.baseline.data_WithOutlier{i}.N); % syls in sequence.
%     
%     for j=1:length(numsyls_max);
%         
%         % === FOR THIS SYLLALBE, plot scatter against pre and post motifs,
%         % for all other syls
%         
%         ffval_thissyl=[];
%         ffval_othersyl=[];
%         
%         
%     
%     lt_figure; hold on;
%     regexpr_string=Params.RegExpr.expressions{i};
%     
%     num_subclasses=length(Params.RegExpr.subexpressions{i});
%     
%     % For each subclass:
%     for ii=1:num_subclasses;
%        
%         subclass=Params.RegExpr.subexpressions{i}{ii};
% 
% 
%         
        
%%
    
    
    
    
    
    
    
%     
%     NumRepeatClasses=length(AllDays_Repeats.RegExpr_string{i}.RepClass);
%     for ii=1:NumRepeatClasses;
%         
%         % === One plot across days for each repeat class
%         repeatclass=AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.string;
%         
%         lt_figure; hold on;
%         title([regexpr_string '; ' repeatclass])
%         
%         for j=1:NumDays;
%             
%             if isempty(AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.FFvals{j})
%                 continue;
%             end
%             
%             Tvals=AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.Tvals_ThisRepClass{j};
%             FFvals=AllDays_Repeats.RegExpr_string{i}.RepClass{ii}.FFvals{j};
%             
%             % convert Tvals to days
%             eventtimes=lt_convert_EventTimes_to_RelTimes(FirstDay, Tvals);
%             Tvals=eventtimes.FinalValue;
%             
%             plot(Tvals, FFvals, 'o');
%             
%             % == PLOT day mean and SEM
%             FFmean=mean(FFvals);
%             FFsem=lt_sem(FFvals);
%             
%             % remove nan
%             FFmean(isnan(FFmean))=[];
%             FFsem(isnan(FFsem))=[];
%             
%             errorbar(ceil(Tvals(1))*ones(length(FFmean),1), FFmean, FFsem, 'o');
%             
%             
%             % ==== PLOT bar graph, to look at change relative to syl repeat
%             % number
%             
%             
%         end
%         
%         
%         
%         SylLabels=AllDays_Repeats.RegExpr_string{i}.RepClass{ii};
%         
%         
%     end
% end


%% SAVE
if saveON==1;
    timestampSv=lt_get_timestamp(0);
    cd(Params.SeqFilter.savedir);
    
    % 
    save('Params','Params');
    save('AllDays_RegExpr','AllDays_RegExpr');
    
    % write a text file that tells you when files were made
    fid1=fopen(['DONE_PlotRegExpr_' timestampSv '.txt'],'w');
    fclose(fid1);
    
    % save figures
    try
        cd FIGURES/PlotRegExpr
    catch err
        mkdir FIGURES/PlotRegExpr;
        cd FIGURES/PlotRegExpr;
    end
    
    lt_save_all_figs
    
    cd ../../
end





end




%% VARIOUS SUBFUNCTIONS
function Fn_AnnotateWNLines(plotWNdays,ylim)

global WNTimeOnInd
global WNTimeOffInd
global DaysToMarkInds

if plotWNdays==1;
    line([WNTimeOnInd-0.5 WNTimeOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
    line([WNTimeOffInd+0.5 WNTimeOffInd+0.5],ylim,'LineStyle','--','Color','r')
    
    for i=1:length(DaysToMarkInds);
        line([DaysToMarkInds{i} DaysToMarkInds{i}],ylim,'LineStyle','--','Color','k');
    end
end
end
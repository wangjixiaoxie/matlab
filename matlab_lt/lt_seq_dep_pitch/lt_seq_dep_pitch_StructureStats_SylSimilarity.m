%% LT 2/5/15 - given feature vectors calculated with lt_seq_dep_pitch_StructureStats, calculates and plots syl similarity scores
function [Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats_SylSimilarity(Params, AllDays_StructStatsStruct)

AllSylFields=fieldnames(AllDays_StructStatsStruct.IndivSyls);
BaselineDays=Params.SeqFilter.BaselineDays;


%% LOOK AT CORRELATIONS IN PITCH - rend to rend

% VISUALIZE BASELINE - smoothed
figure; hold on;
title('Pitch vs. time during baseline, for all syls');
smbin=10; % rends to smooth over

syllist=Params.SeqFilter.SylLists.SylsSame;
PlotColors=lt_make_plot_colors(length(syllist),0);

for i=1:length(syllist);
    syl=syllist{i};
    
    % extract baseline Inds
    bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,BaselineDays); % inds of rends that are from baseline
    
%     times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum(bInds); % times for all rends
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(:,2);
    vals_sm=smooth(vals,smbin);
   
    
      hplot(i)=plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});

    
%     if sum(strcmp(syl,Params.SeqFilter.SylLists.TargetSyls))==1; % color based on relationship to targ syl.
%         
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %          plot(times.FinalValue,vals_sm+200,'-k');
%        
%     elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsSame))==1;
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %        plot(times.FinalValue,vals_sm+200,'-b');
% 
%     elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsDifferent))==1;
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %        plot(times.FinalValue,vals_sm+200,'-r');
%         
%     else
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %        plot(times.FinalValue,vals_sm+200,'-y');
%         disp(['Problem - ' syl ' is not targ, similar, or diff. what is it?']);
%     end
end

% plot targ
syl=Params.SeqFilter.SylLists.TargetSyls{1};
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(:,2);
    vals_sm=smooth(vals,smbin);
      hplot(i+1)=plot(times.FinalValue,vals_sm+200,'-','Color','k');

legend(hplot,syllist)

figure
for i=1:length(AllSylFields);
    htmp(i)=plot(1,1,'Color',PlotColors{i});
end
legend(htmp,AllSylFields)

%% LOOK AT FLUCTUATION OVER TIME - 1ST pass, deviation from day mean.


% Z-SCORE all data versus day mean/std
% all renditions, get z-score of feature vector relative to the
% mean/std feature vector of that day.
syllist=AllSylFields;
NumDays=datenum(Params.SeqFilter.LastDay)-datenum(Params.SeqFilter.FirstDay)+1;

for i=1:length(syllist);
    syl=syllist{i};
    
    % for each day, subtract mean of that day
    for ii=1:NumDays;
        
        inds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,ii);
        FVdayvals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(inds,:);
        FVdaymean=mean(FVdayvals);
        FVdaystd=std(FVdayvals);
        
        FVdayzscores=(FVdayvals-repmat(FVdaymean,size(FVdayvals,1),1))./repmat(FVdaystd,size(FVdayvals,1),1);
        
        AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(inds,:)=FVdayzscores;
        
    end
end


% PLOT over all days (PITCH)
figure; hold on; 

syllist=Params.SeqFilter.SylLists.SylsSame;
PlotColors=lt_make_plot_colors(length(syllist),0);
smbin=20;


for i=1:length(syllist);
    syl=syllist{i};
    

times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
    vals=mean(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean,2); % take mean z-score of all vals
%     vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(:,2);
    vals_sm=smooth(vals,smbin);
   
    
      hplot(i)=plot(times.FinalValue,vals_sm+0.2*i,'-','Color',PlotColors{i});
%       hplot(i)=plot(times.FinalValue,vals_sm,'.','Color',PlotColors{i});

    
%     if sum(strcmp(syl,Params.SeqFilter.SylLists.TargetSyls))==1; % color based on relationship to targ syl.
%         
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %          plot(times.FinalValue,vals_sm+200,'-k');
%        
%     elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsSame))==1;
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %        plot(times.FinalValue,vals_sm+200,'-b');
% 
%     elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsDifferent))==1;
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %        plot(times.FinalValue,vals_sm+200,'-r');
%         
%     else
%         plot(times.FinalValue,vals_sm+200,'-','Color',PlotColors{i});
% %        plot(times.FinalValue,vals_sm+200,'-y');
%         disp(['Problem - ' syl ' is not targ, similar, or diff. what is it?']);
%     end
end

% plot targ
syl=Params.SeqFilter.SylLists.TargetSyls{1};
    times=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum; % times for all rends
    times=lt_convert_EventTimes_to_RelTimes(Params.SeqFilter.FirstDay,times); % output structure with times
%     vals=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean(:,2);
    vals=mean(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect_ZscoresRelDayMean,2); % take mean over all feature z-scores
    
    vals_sm=smooth(vals,smbin);
      hplot(i+1)=plot(times.FinalValue,vals_sm,'-','Color','k');
%       hplot(i+1)=plot(times.FinalValue,vals_sm,'.','Color','k');

legend(hplot,syllist)






%% GET MEAN BASELINE VECTOR

for i=1:length(AllSylFields);
    syl=AllSylFields{i};
    
    bInds=ismember(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd,BaselineDays); % inds of rends that are from baseline
    fvec=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect(bInds,:); % these are baseline, all renditions
    
    % mean and std
    AllFeatVecMean(i,:)=mean(fvec,1);
    AllFeatVecSD(i,:)=std(fvec);
    AllFeatVecN(i)=size(fvec,1);
    
end



%% FIRST, FOR EACH FEATURE, look at values for all syls
% ONLY LOOK AT BASELINE


% PLOT - one plot for each feature vector, compare syls;
for i=1:length(Params.StructureStats.FeatureLegend); % for each feature
    feature=Params.StructureStats.FeatureLegend{i};
    
    figure; hold on;
    title(['Mean (SD) of ' feature ' (during baseline)']);
    
    for ii=1:length(AllSylFields);
        syl=AllSylFields{ii};
        
        if sum(strcmp(syl,Params.SeqFilter.SylLists.TargetSyls))==1; % color based on relationship to targ syl.
            
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'ok','MarkerSize',9);
        elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsSame))==1;
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'og','MarkerSize',9);
        elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsDifferent))==1;
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'or','MarkerSize',9);
        else
            errorbar(ii,AllFeatVecMean(ii,i),AllFeatVecSD(ii,i),'ob','MarkerSize',9);
            disp(['Problem - ' syl ' is not targ, similar, or diff. what is it?']);
        end
        
    end
    
    set(gca,'XTick',1:length(AllSylFields));
    
    set(gca,'XTickLabel',AllSylFields);
    
end






%% SAVE








function [SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_SummaryPlot(SeqDepPitch_AcrossBirds, DayBinSize)
%% LT 5/13/15 - use as a subfunction of lt_seq_dep_pitch_ACROSSBIRDS
% Plot each experiment/bird separately.

% Also get's some stats to put into the data structure, which coudl come in
% handy later


%% Params

NumBirds_OneTarg=length(SeqDepPitch_AcrossBirds.birds);


%% RAW PLOTS for each experient (contains consolidation window as well)


% First, plot raw data for all experiments, annotating consolidation
% Second, save indices for consolid start and end.

YTOT={};
for i=1:NumBirds_OneTarg;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylList=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions);
        
        plotcols=lt_make_plot_colors(length(SylList),0,0);
        lt_figure; hold on;
        
        
        hplot=[];
        YTOT{i,ii}.beginning=[];
        YTOT{i,ii}.end=[];
        for iii=1:length(SylList);
            syl=SylList{iii};
            
            % == SUBPLOT 1 - line plots
            lt_subplot(2,1,1); hold on;

            % == Put lines for consolidation start and end
            ConsolStartDate=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartDate;
            ConsolEndDate=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndDate;
            FirstDate=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
            
            % convert dates to indices
            ConsolStartInd=datenum(ConsolStartDate,'ddmmmyyyy')-datenum(FirstDate,'ddmmmyyyy')+1;
            ConsolEndInd=datenum(ConsolEndDate,'ddmmmyyyy')-datenum(FirstDate,'ddmmmyyyy')+1;
            
            % --- get values to plot
            FF_DayVals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase;
            
            plot(1:length(FF_DayVals),FF_DayVals,'-o','Color',plotcols{iii},'LineWidth',2);
            
            
            % === Plot dots for pitch at beginning and end of consolidation
            % get mean for beginning and end of consolid 
           
            FF_consol_start=nanmean(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase(ConsolStartInd:ConsolStartInd+DayBinSize-1));
            FF_consol_end=nanmean(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase(ConsolEndInd-DayBinSize+1:ConsolEndInd));
            
            % === OUTPUT Structures
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_FF_minusBase=FF_consol_start;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolEnd_FF_minusBase=FF_consol_end;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndInd=ConsolEndInd;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd=ConsolStartInd;
            
            % === save to get global means
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                % then is target, save to target field
                YTOT{i,ii}.target_beginning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolStart_FF_minusBase;
                YTOT{i,ii}.target_end=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).ConsolEnd_FF_minusBase;
            else
                % then is not target, save to general field;
                YTOT{i,ii}.beginning=[YTOT{i,ii}.beginning, FF_consol_start];
                YTOT{i,ii}.end=[YTOT{i,ii}.end, FF_consol_end];
            end
            % =====
            
            
            % plot (one dot for each syl)
            hplot(iii)=plot(ConsolStartInd+DayBinSize/2, FF_consol_start, 's', 'Color','k','MarkerSize',8,'MarkerFaceColor',plotcols{iii});
            plot(ConsolEndInd-DayBinSize/2, FF_consol_end ,'s', 'Color','k','MarkerSize',8,'MarkerFaceColor',plotcols{iii});
            
            % plot lines marking epoch start and end
            line([ConsolStartInd-0.5 ConsolStartInd-0.5],ylim);
            line([ConsolEndInd+0.5 ConsolEndInd+0.5],ylim);
            
            lt_plot_zeroline;            

        end
        
        legend(hplot, SylList)
        
        % ----------- Put lines for WN start and end, and baseline
        WNOnInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        WnOffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        baseline_last_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays(end);
        line([WNOnInd-0.5 WNOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
        line([WnOffInd+0.5 WnOffInd+0.5],ylim,'LineStyle','--','Color','r')
        line([baseline_last_day+0.4 baseline_last_day+0.4],ylim,'LineStyle','--','Color','g')
        
        % == SUBPLOT 2 - just dots for beginning and end (not absolute values);
        lt_subplot(2,3,4); hold on;
        title('learning (start and end of consolidation)');
        
        X=[1.1, 2.1];
        
        % plot nontargets
        for j=1:length(YTOT{i,ii}.beginning);
            Y1=YTOT{i,ii}.beginning(j);
            Y2=YTOT{i,ii}.end(j);
            
            lt_plot(X, [Y1 Y2],{'LineStyle','-','Color',plotcols{j}});
            
            xlim([-0.5 3.5]);
        end
        
        % plot the target differently
        Y1=YTOT{i,ii}.target_beginning;
        Y2=YTOT{i,ii}.target_end;
        
        plot(X, [Y1 Y2],'-o','Color','k');
        
        % -- add means
        Ybeg_mean=mean(YTOT{i,ii}.beginning);
        Ybeg_sem=lt_sem(YTOT{i,ii}.beginning);
        
        Yend_mean=mean(YTOT{i,ii}.end);
        Yend_sem=lt_sem(YTOT{i,ii}.end);
        
        errorbar(X+0.2, [Ybeg_mean Yend_mean], [Ybeg_sem Yend_sem],'-s', 'MarkerFaceColor','k','MarkerSize',9);
        
        lt_plot_zeroline;
        
        
        % === SUBPLOT 3 -  absolute values
        lt_subplot(2,3,5); hold on;
        title('absolute value of learning (consolidaton)');
        
        X=[1.1, 2.1];
        
        % plot nontargets
        for j=1:length(YTOT{i,ii}.beginning);
            Y1=abs(YTOT{i,ii}.beginning(j));
            Y2=abs(YTOT{i,ii}.end(j));
            
            lt_plot(X, [Y1 Y2],{'LineStyle','-','Color',plotcols{j}});
            
            xlim([-0.5 3.5]);
        end
        
        % plot the target differently
        Y1=abs(YTOT{i,ii}.target_beginning);
        Y2=abs(YTOT{i,ii}.target_end);
        
        plot(X, [Y1 Y2],'-o','Color','k');
        
        % -- add means
        Ybeg_mean=mean(abs(YTOT{i,ii}.beginning));
        Ybeg_sem=lt_sem(abs(YTOT{i,ii}.beginning));
        
        Yend_mean=mean(abs(YTOT{i,ii}.end));
        Yend_sem=lt_sem(abs(YTOT{i,ii}.end));
        
        errorbar(X+0.2, [Ybeg_mean Yend_mean], [Ybeg_sem Yend_sem],'-s', 'MarkerFaceColor','k','MarkerSize',9);


        % === SUBPLOT 4 - maginute of learning (recalculated at start and
        % end)
        subplot(2,3,6); hold on;
        title('learning magnitude (relative to target)');
        
        X=[1.1, 2.1];
        
        % plot nontargets
        for j=1:length(YTOT{i,ii}.beginning);
            Y1=YTOT{i,ii}.beginning(j)/YTOT{i,ii}.target_beginning;
            Y2=YTOT{i,ii}.end(j)/YTOT{i,ii}.target_end;
            
            lt_plot(X, [Y1 Y2],{'LineStyle','-','Color',plotcols{j}});
            
            xlim([-0.5 3.5]);
        end
        
        
        % -- add means
        
        Ybeg_mean=mean(YTOT{i,ii}.beginning./YTOT{i,ii}.target_beginning);
        Ybeg_sem=lt_sem(YTOT{i,ii}.beginning./YTOT{i,ii}.target_beginning);
        
        Yend_mean=mean(YTOT{i,ii}.end./YTOT{i,ii}.target_end);
        Yend_sem=lt_sem(YTOT{i,ii}.end./YTOT{i,ii}.target_end);
        
        errorbar(X+0.2, [Ybeg_mean Yend_mean], [Ybeg_sem Yend_sem],'-s', 'MarkerFaceColor','k','MarkerSize',9);

        
        % === Annotate
        lt_subtitle([birdname '; experiment: ' exptname]);
        
    end
end


% ==== PLOT ACROSS ALL BIRDS
lt_figure; hold on;
plotcols=lt_make_plot_colors(size(YTOT,1),0,0);
title('During consolidation, All syls, all birds');
xlabel('Start');
ylabel('End');

Xtot2=[];
Ytot2=[];
GroupIndicator=[];

for i=1:size(YTOT,1); % each bird
    for ii=1:length(YTOT(i,:)); % bird's experiemnts
        
        % -- scatterplot of all raw data (learning norm to target)
        if ~isempty(YTOT{i,ii})
        X= YTOT{i,ii}.beginning./YTOT{i,ii}.target_beginning;
        Y= YTOT{i,ii}.end./YTOT{i,ii}.target_end;
        
        plot(X,Y,'o','Color',plotcols{i});
        
        
        % -- collect all values across all birds
        if isnan(X);
            keyboard
        end
        Xtot2=[Xtot2 X];
        Ytot2=[Ytot2 Y];
        
        
        
        end
    end
end

xlim([-0.8, 0.8]);
ylim([-0.8 0.8]);
      
line(xlim,[0 0])
line([0 0], ylim);
line(xlim,ylim,'Color','k')

xlabel('Learning (beginning)');
ylabel('Learning (end)');


% ==== PLOT DOTS
lt_figure; hold on;
title('Consolidation (start and end)');
ylabel('Learning (rel target)');

for i=1:length(Xtot2);
    plot([1 2], [Xtot2(i), Ytot2(i)],'-o'); 
end
xlim([-0.5 3.5]);

% plot mean
Y=[mean(Xtot2), mean(Ytot2)];
Ysem=[lt_sem(Xtot2) lt_sem(Ytot2)];
errorbar([1.1 2.1], Y, Ysem, 'sk','MarkerFaceColor','k','MarkerSize',8);

lt_plot_zeroline;

% -- ttest
% [~, p]=ttest(Xtot2,Ytot2);
[p,~]=signrank(Xtot2,Ytot2);
text(1.5,0.8,['sign rank p=' num2str(p)],'FontSize',13,'FontWeight','bold');


% ==== PLOT DOTS (absolute)
lt_figure; hold on;
title('Absolute learning');


Xtot2_abs=abs(Xtot2);
Ytot2_abs=abs(Ytot2);


for i=1:length(Xtot2);
    plot([1 2], [Xtot2_abs(i), Ytot2_abs(i)],'-o'); 
end
xlim([-0.5 3.5]);

% plot mean
Y=[mean(Xtot2_abs), mean(Ytot2_abs)];
Ysem=[lt_sem(Xtot2_abs) lt_sem(Ytot2_abs)];
errorbar([1.1 2.1], Y, Ysem, 'sk','MarkerFaceColor','k','MarkerSize',8);

lt_plot_zeroline;

% -- ttest
[~, p]=ttest(Xtot2_abs,Ytot2_abs);
[p,~]=signrank(Xtot2_abs,Ytot2_abs);
text(1.5,0.8,['sign rank p=' num2str(p)],'FontSize',13,'FontWeight','bold');



% === PLOT CHANGE (over consolidation) relative to learning
% -- Get data in form of difference (i.e. end minus beginning)
Y_EndMinusStart=Ytot2-Xtot2;

% plot
lt_figure; hold on;

% regression
[~,~,~,~,~,SummaryStats]=lt_regress(Y_EndMinusStart,Xtot2,1);

title('Change over consolidation vs. starting learning');
ylabel('End minus start (generalization)');
xlabel('Learning at start of consolidation');

line(xlim,[0 0]);
line([0 0],ylim);







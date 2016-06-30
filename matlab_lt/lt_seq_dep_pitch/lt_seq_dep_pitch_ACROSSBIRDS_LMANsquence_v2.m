function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANsequence_v2(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 12/11/15 - not using RegExp code. instead, enter variable branch point and will go thru and find them in raw dat
% NOTE: ASSUMES THAT ALL SONGS TO ANALYZE HAS "B" LABELED 

count=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


assert(mod(length(PARAMS.sequence_v2.birdname_exptname_Motifs), 3)==0, 'Error, did not enter params correctly (mod3 should be 0)');
% KeySylToExtractSongLabels='Ab'; % a key syl that should be present in all songs that want to quantify
ConditionFieldList={'PBS','MUSC'};

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


if input('WARNING - code takes long time. still want to run? (1 or 0)')==0
    asdfassv
end

%%

for i=1:NumBirds;
    
    NumExpt=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:NumExpt
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % == Search input to see if this expt shoud be plotted
        
        ind1=find(strcmp(PARAMS.sequence_v2.birdname_exptname_Motifs, birdname));
        ind2=find(strcmp(PARAMS.sequence_v2.birdname_exptname_Motifs, exptname));
        
        if isempty(ind1) || isempty(ind2);
            continue
        else
            if length(intersect(ind1, ind2-1))==0;
                continue
            elseif length(intersect(ind1, ind2-1))>1;
                disp('PROBLEM... Params not defined properly');
                fasfasdceafa
            end
        end
        
        
        % === RUN ANALYSIS AND PLOT FOR THIS BIRRD EXPT
        ind=intersect(ind1, ind2-1)+2;
        MotifsToCount= PARAMS.sequence_v2.birdname_exptname_Motifs{ind};
        disp(MotifsToCount);
        
        
        
        %% LOAD
        
        DirOfRawDat=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/AllDays_RawDatStruct.mat'];
        DirOfRawParams=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir '/Params.mat'];
        
        load(DirOfRawDat);
        load(DirOfRawParams);
        
        
        
        %% === For each day, get probability for each transition
        % DaysToCheck=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
        targetsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targetsyl).FFvals);
        DaysToCheck=1:NumDays;
        
        DATSTRUCT=struct;
        
        
        for day=DaysToCheck
            
            % ================================ PBS - figure out what are raw datapoint to use
            if ~isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targetsyl).Tvals_WithinTimeWindow{day})
                
                % Extract all potential tvals of songs in this time window
                SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                TvalsWithinTimeWindow=[];
                for ss=1:length(SylsUnique);
                    syl=SylsUnique{ss};
                    
                    Tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                    TvalsWithinTimeWindow=[TvalsWithinTimeWindow Tvals];
                end
                TvalsWithinTimeWindow=unique(TvalsWithinTimeWindow);
                %
                %         TvalsWithinTimeWindow=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(KeySylToExtractSongLabels).Tvals_WithinTimeWindow{day};
                TimeAllRawDat=cell2mat(AllDays_RawDatStruct{day}.data_WithOutlier.b(:,6));
                
                % -- find overlap of raw times and good times
                [~, indA, indB]=intersect(TvalsWithinTimeWindow, TimeAllRawDat);
                
                % Make sure there are exactly correct number of songs extracted
                assert(length(indB)==length(unique(TvalsWithinTimeWindow)), 'PROBLEM, number of songs extracted does not match expected number');
                
                % EXTRACT song labels - only one per song
                SongLabels_unique=AllDays_RawDatStruct{day}.data_WithOutlier.b(indB, 7);
                
            else
                SongLabels_unique=[];
            end
            
            % Count total transitiosn
            TotalAcrossMotifs=0;
            
            % Count number of transition of desired type
            for j=1:length(MotifsToCount);
                Motif=MotifsToCount{j};
                
                
                NumOccurances_all=0;
                NumSongs=length(SongLabels_unique);
                
                for k=1:length(SongLabels_unique);
                    SongLabel=SongLabels_unique{k};
                    
                    NumOccurances=findstr(SongLabel, Motif);
                    
                    NumOccurances_all=NumOccurances_all+length(NumOccurances);
                    
                end
                
                % ==== OUTPUT
                DATSTRUCT.data.PBS.NumOccurances(day).(Motif)=NumOccurances_all;
                DATSTRUCT.data.PBS.NumSongs(day)=NumSongs;
                
                TotalAcrossMotifs=TotalAcrossMotifs+NumOccurances_all;
            end
            
            DATSTRUCT.data.PBS.TotalAcrossMotifs(day)=TotalAcrossMotifs;
            
            % ================================ MUSC - figure out what are raw datapoint to use
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning, 'DataMatrix_MUSC')
                if ~isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targetsyl).Tvals_WithinTimeWindow{day})
                    
                    % get tvals of all songs in this condition window
                    SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    TvalsWithinTimeWindow=[];
                    for ss=1:length(SylsUnique);
                        syl=SylsUnique{ss};
                        
                        Tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                        TvalsWithinTimeWindow=[TvalsWithinTimeWindow Tvals];
                    end
                    TvalsWithinTimeWindow=unique(TvalsWithinTimeWindow);
                    
                    
                    %             TvalsWithinTimeWindow=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(KeySylToExtractSongLabels).Tvals_WithinTimeWindow{day};
                    TimeAllRawDat=cell2mat(AllDays_RawDatStruct{day}.data_WithOutlier.b(:,6));
                    
                    % -- find overlap of raw times and good times
                    [~, indA, indB]=intersect(TvalsWithinTimeWindow, TimeAllRawDat);
                    
                    % EXTRACT song labels - only one per song
                    SongLabels_unique=AllDays_RawDatStruct{day}.data_WithOutlier.b(indB, 7);
                else
                    
                    SongLabels_unique=[];
                end
                
                % Count total transitions
                TotalAcrossMotifs=0;
                
                % Count number of transition of desired type
                for j=1:length(MotifsToCount);
                    Motif=MotifsToCount{j};
                    
                    
                    NumOccurances_all=0;
                    NumSongs=length(SongLabels_unique);
                    
                    for k=1:length(SongLabels_unique);
                        SongLabel=SongLabels_unique{k};
                        
                        NumOccurances=findstr(SongLabel, Motif);
                        
                        NumOccurances_all=NumOccurances_all+length(NumOccurances);
                        
                    end
                    
                    % ==== OUTPUT
                    DATSTRUCT.data.MUSC.NumOccurances(day).(Motif)=NumOccurances_all;
                    DATSTRUCT.data.MUSC.NumSongs(day)=NumSongs;
                    TotalAcrossMotifs=TotalAcrossMotifs+NumOccurances_all;
                    
                end
                DATSTRUCT.data.MUSC.TotalAcrossMotifs(day)=TotalAcrossMotifs;
                
            end
        end
        
        
        %% === PLOT, for each day, plot motif probability
        
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);

        plotcols=lt_make_plot_colors(length(MotifsToCount),0,0);
        
        
        for l=1:length(ConditionFieldList);
            
            ConditionField=ConditionFieldList{l};
            
            
            
            for k=1:length(MotifsToCount);
                Motif=MotifsToCount{k};
                
                if DATSTRUCT.data.(ConditionField).NumSongs(j)<0
                    Yprob=[Yprob nan];
                    continue
                end
                
                % ==== calc probability for this motif\
                DayInds=find(DATSTRUCT.data.(ConditionField).TotalAcrossMotifs~=0);
                
                Yprob=[DATSTRUCT.data.(ConditionField).NumOccurances(DayInds).(Motif)]./DATSTRUCT.data.(ConditionField).TotalAcrossMotifs(DayInds);
                
                % === PLOT
                if strcmp(ConditionField, 'MUSC');
                    lt_plot(DayInds, Yprob, {'Color', plotcols{k},  'Marker','s', 'MarkerFaceColor','none'});
                else
                    lt_plot(DayInds, Yprob, {'Color', plotcols{k}, 'LineStyle', '-'});
                end
            end
            
            % --- PLOT LINES FOR EVENTS
            % ============== annotate with lines for stuff
            Ylim=ylim;
            Yrange=Ylim(2)-Ylim(1);
            Xlim=xlim;
            
            % WN on and off
            WNonday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            line([WNonday-0.5 WNonday-0.5], ylim, 'Color' ,'r');
            WNoffday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            line([WNoffday+0.5 WNoffday+0.5],ylim,'Color','r');
            lt_plot_text(WNonday, Ylim(1), 'WNon/off','r')
            
            % day I called consol start
            try
                day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
                day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndInd;
                
                line([day1-0.5 day1-0.5], ylim);
                line([day2+0.5 day2+0.5], ylim);
                
                lt_plot_text(day1, Ylim(1), 'consol','b')
            catch err
            end
            
            % multidir
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
                date1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
                date2_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
                
                line([date1_ind+0.4 date1_ind+0.4], ylim, 'Color', 'k');
                line([date2_ind+0.4 date2_ind+0.4], ylim, 'Color', 'k');
                lt_plot_text(date1_ind, Ylim(1), 'multidir','k')
            end
            
            %         % LMAN inactivation days
            %         if  SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            %             targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            %             muscdays_inds=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow));
            %
            %             % plot
            %             plot(muscdays_inds, Ylim(1)+Yrange/6, '^k', 'MarkerSize', 8);
            %             plot(muscdays_inds, Ylim(2)-Yrange/6, 'vk', 'MarkerSize', 8);
            %
            %             %             for k=1:length(muscdays_inds);
            %             % %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'Color',[0.7 0.7 0.7]);
            %             %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'LineStyle','--','Color',[0.3 0.3 0.3]);
            %             %             end
            %
            %             lt_plot_text(1, Ylim(1)+Yrange/8, 'arrowhead: MUSC', 'k');
            %
            %
            %         end
            
            %         % annotate target learning direction (as detected automatically)
            %         targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            %
            %         lt_plot_text(2, Ylim(2)-Yrange/8, ['targ learn dir: ' num2str(targ_learn_dir)], 'r');
            
            % ---- what was target?
            target=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            lt_plot_text(0,Ylim(2)-Yrange/8, ['target: ' target]);
            
            try
                multi_targs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
                
                lt_plot_text(0,0.8, ['multidir targs: ' multi_targs]);
            catch err
            end
            
            
        end
        
        legend(MotifsToCount);
        
        
    end
end


function lt_seq_dep_pitch_ACROSSBIRDS_LMANRepeats(SeqDepPitch_AcrossBirds, Params, epoch_field)
%% LT 11/25/15 -
% epoch_field='days_consolid_early';

plot_text_names=0;
DoNotKeepNormDataIfTargLearnLessThanThresh=1;
learn_thresh=0.5; % zscore



%% EXTRACT ONLY EXPERIMENTS WITH REPEATS
SeqDepPitch_AcrossBirds = lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, 'RepeatExptsOnly');

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% EXTRACT EXPERIMENTS WITH ONLY 
%% COLLECT DATA: output is  {targlocation}{repeatlength}[N x SylInRepeat]

RepeatLearning_ByTargLocation=cell(1,4); % {targlocation}{repeatlength}[N x SylInRepeat]
RepeatMPbias_ByTargLocation=cell(1,4);

% 
% RepeatFFsem_ByTargLocation=cell(1,4); % {targlocation}{repeatlength}[N x SylInRepeat]
% 
% RepeatFFmeanNORM_ByTargLocation=cell(1,4);
% RepeatFFsemNORM_ByTargLocation=cell(1,4);
% 
BirdExptTargNames=cell(1,4);

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % ==== check where target hit
        targpos=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DataRepeats.PosOfTargInRepeat;
        
        % ==== ffvals for all pos in repeat
        SylsInRepeat=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DataRepeats.SylsInRepeat;
        replength=length(SylsInRepeat);
        
        learning_all=[];
        MPbias_all=[];

        for j=1:length(SylsInRepeat);
            syl=SylsInRepeat{j};
            
            %             days=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DayInds;
            if any(strcmp(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique, syl))
                % then this syl has data
                if isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).mean_Zscore);
                    replength=j-1;
                    break % i.e. end of repeat
                end
                
                % ===== GET STATS
%                 zscore=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).mean_Zscore;
%                 zscore_all=[zscore_all zscore];
%                 
%                 zscore_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).sem_Zscore;
%                 zscore_sem_all=[zscore_sem_all zscore_sem];
                
                % ------ LMAN STATS
                learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).meanFF_pbs;
                learning_all=[learning_all learning];
                
                MPbias=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).meanFF_musc;
                MPbias_all=[MPbias_all MPbias];
               
                
%                 learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).semFF_pbs;
%                 MPbias_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.(epoch_field).(syl).semFF_musc;
                
            else
                replength=j-1;
                break
            end
        end
        
        % ==== throw out if replength is <2
        if replength<2
            disp(['THROWING OUT: ' birdname '-' exptname ' (replength<2)']);
            continue
        end
           
        
        % ==== flip sign if learning is neg at targ
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
            learning_all=-learning_all;
            
            MPbias_all=-MPbias_all;
            
        end
        
%             % normalize all syls to targ syl
%             % put nans if targ syl learned less than 0.5z (optional)
%             zscore_all_norm=zscore_all./zscore_all(targpos);
%             zscore_sem_all_norm=zscore_sem_all./zscore_all(targpos);
%             
%             if DoNotKeepNormDataIfTargLearnLessThanThresh==1;
%                 if zscore_all(targpos)<learn_thresh;
%                     zscore_all_norm=nan(1, length(zscore_all));
%                     zscore_sem_all_norm=nan(1, length(zscore_all));
%                 end
%             end
% 
%             
% %         if targpos==2 && replength==4;
% %             keyboard
% %         end
%         
        % ==== put into output
        if length(RepeatLearning_ByTargLocation{targpos})<replength;
            RepeatLearning_ByTargLocation{targpos}{replength}=learning_all;
        else
            RepeatLearning_ByTargLocation{targpos}{replength}=[RepeatLearning_ByTargLocation{targpos}{replength}; learning_all];
        end
        
        
                % ==== put into output
        if length(RepeatMPbias_ByTargLocation{targpos})<replength;
            RepeatMPbias_ByTargLocation{targpos}{replength}=MPbias_all;
        else
            RepeatMPbias_ByTargLocation{targpos}{replength}=[RepeatMPbias_ByTargLocation{targpos}{replength}; MPbias_all];
        end
        
        
        %         % ==== put into output
%         if length(RepeatFFsem_ByTargLocation{targpos})<replength;
%             RepeatFFsem_ByTargLocation{targpos}{replength}=zscore_sem_all;
%         else
%             RepeatFFsem_ByTargLocation{targpos}{replength}=[RepeatFFsem_ByTargLocation{targpos}{replength}; zscore_sem_all];
%         end
%         
%         % ==== put into output
%         if length(RepeatFFmeanNORM_ByTargLocation{targpos})<replength;
%             RepeatFFmeanNORM_ByTargLocation{targpos}{replength}=zscore_all_norm;
%         else
%             RepeatFFmeanNORM_ByTargLocation{targpos}{replength}=[RepeatFFmeanNORM_ByTargLocation{targpos}{replength}; zscore_all_norm];
%         end
%       
%         % ==== put into output
%         if length(RepeatFFsemNORM_ByTargLocation{targpos})<replength;
%             RepeatFFsemNORM_ByTargLocation{targpos}{replength}=zscore_sem_all_norm;
%         else
%             RepeatFFsemNORM_ByTargLocation{targpos}{replength}=[RepeatFFsemNORM_ByTargLocation{targpos}{replength}; zscore_sem_all_norm];
%         end
%       
%         
%         
%         
        % =====
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        if length(BirdExptTargNames{targpos})<replength;
            stringtoadd=[birdname '_' exptname '_' targsyl];
            BirdExptTargNames{targpos}{replength}={stringtoadd};
        else
            BirdExptTargNames{targpos}{replength}=[BirdExptTargNames{targpos}{replength} {[birdname '_' exptname '_' targsyl]}];
        end
        
    end
end


%% PLOT LEARNING AND MP BIAS BY [targ location] and [rep class] [NON NORMALIZED]


count=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


NumTargLocations=length(RepeatLearning_ByTargLocation);

hsplots=[];

cc=1;
for i=1:NumTargLocations;
    
    if isempty(RepeatLearning_ByTargLocation{i})
        continue;
    end
    
    NumRepClassesActual=sum(~cellfun(@isempty, RepeatLearning_ByTargLocation{i}));
    
    %     for ii=1:NumRepClasses;
    
    NumClassesOfRepeatLengths=length(RepeatLearning_ByTargLocation{i});
    
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatLearning_ByTargLocation{i}{ii})
            continue
        end
        
%         hsplot=lt_subplot(NumTargLocations, NumRepClassesActual, cc);
        [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);

        title(['targ loc: ' num2str(i), '; rep class: ' num2str(ii)]);
        
        
        hold on;
        hsplots=[hsplots hsplot];
        
        % ====== PLOT INDIVIUDAL EXPERIMENTS
        x=1:ii;
        
        for k=1:size(RepeatLearning_ByTargLocation{i}{ii},1);
            %             lt_plot(x, RepeatLearning_ByTargLocation{i}{ii}(k,:), {'Errors', RepeatFFsem_ByTargLocation{i}{ii}(k,:), 'LineStyle','-', 'MarkerFaceColor','none'});
            
            % ------ LEARNING
            Y=[RepeatLearning_ByTargLocation{i}{ii}(k,:); RepeatMPbias_ByTargLocation{i}{ii}(k,:)]; % [learning; MPbias];
            
            for kk=1:size(Y, 2);
                plot([kk-0.1 kk+0.1], Y(:,kk), '-k');
                
                plot([kk-0.1], Y(1,kk), 'ok'); % learning
                plot([kk+0.1], Y(2,kk), 'or'); % Mp bias
                
                try
                    plot([kk-0.1 kk+0.9], [Y(1,kk) Y(1, kk+1)], '--'); % connectign line to next rep in repeat
                catch err
                    % reachedc end of repeat
                end
                
            end
            
            if plot_text_names==1;
                string_to_plot=BirdExptTargNames{i}{ii}{k};
                
                text(x(end), RepeatLearning_ByTargLocation{i}{ii}(k, end), string_to_plot);
                
            end
            
        end
        
        
        % ===== plot mean for this replength class
        if size(RepeatLearning_ByTargLocation{i}{ii},1)>1;
            
            Ymean_learn=mean(RepeatLearning_ByTargLocation{i}{ii});% learning
            Ysem_learn=lt_sem(RepeatLearning_ByTargLocation{i}{ii});
            
            Ymean_MP=mean(RepeatMPbias_ByTargLocation{i}{ii}); % MP
            Ysem_MP=lt_sem(RepeatMPbias_ByTargLocation{i}{ii});
            
            lt_plot_bar(x-0.15, Ymean_learn, {'Errors', Ysem_learn, 'Color','k', 'BarWidth',0.3})
            lt_plot_bar(x+0.15, Ymean_MP, {'Errors', Ysem_MP, 'Color','r', 'BarWidth',0.3})
            
            
            
        end
        cc=cc+1;
        
        lt_plot_zeroline;
        
        % ==== INDICATE TARGET LOCATION
        plot(i, 0, '^r', 'MarkerSize',10)
        
    end
    %     end
end


% lt_plot_equalyaxis(hsplots);

lt_plot_equalaxes(hsplots,1,1);
% linkaxes(hsplots, 'x','y');





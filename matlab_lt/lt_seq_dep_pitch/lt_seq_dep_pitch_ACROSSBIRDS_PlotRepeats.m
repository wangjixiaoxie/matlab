function lt_seq_dep_pitch_ACROSSBIRDS_PlotRepeats(SeqDepPitch_AcrossBirds, Params)
%% LT 11/25/15 -

plot_text_names=0;
DoNotKeepNormDataIfTargLearnLessThanThresh=1;
learn_thresh=0.5; % zscore

%% EXTRACT ONLY EXPERIMENTS WITH REPEATS
SeqDepPitch_AcrossBirds = lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, 'RepeatExptsOnly');


NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% COLLECT DATA: output is  {targlocation}{repeatlength}[N x SylInRepeat]

TargetLocationList=[1 2 3]; % i.e. filtering out experiemnts where targ was either 1st, 2nd, or .... in repeat

RepeatFFmean_ByTargLocation=cell(1,4); % {targlocation}{repeatlength}[N x SylInRepeat]
RepeatFFsem_ByTargLocation=cell(1,4); % {targlocation}{repeatlength}[N x SylInRepeat]

RepeatFFmeanNORM_ByTargLocation=cell(1,4);
RepeatFFsemNORM_ByTargLocation=cell(1,4);

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
        zscore_all=[];
        zscore_sem_all=[];
        
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
                zscore=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).mean_Zscore;
                zscore_all=[zscore_all zscore];
                
                zscore_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.TargPassLearnThresh.DATA.(syl).sem_Zscore;
                zscore_sem_all=[zscore_sem_all zscore_sem];
                
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
            zscore_all=-zscore_all;
        end
        
                    % normalize all syls to targ syl
            % put nans if targ syl learned less than 0.5z (optional)

            zscore_all_norm=zscore_all./zscore_all(targpos);
            zscore_sem_all_norm=zscore_sem_all./zscore_all(targpos);
            
            if DoNotKeepNormDataIfTargLearnLessThanThresh==1;
                if zscore_all(targpos)<learn_thresh;
                    zscore_all_norm=nan(1, length(zscore_all));
                    zscore_sem_all_norm=nan(1, length(zscore_all));
                end
            end

            
%         if targpos==2 && replength==4;
%             keyboard
%         end
        
        % ==== put into output
        if length(RepeatFFmean_ByTargLocation{targpos})<replength;
            RepeatFFmean_ByTargLocation{targpos}{replength}=zscore_all;
        else
            RepeatFFmean_ByTargLocation{targpos}{replength}=[RepeatFFmean_ByTargLocation{targpos}{replength}; zscore_all];
        end
        
        % ==== put into output
        if length(RepeatFFsem_ByTargLocation{targpos})<replength;
            RepeatFFsem_ByTargLocation{targpos}{replength}=zscore_sem_all;
        else
            RepeatFFsem_ByTargLocation{targpos}{replength}=[RepeatFFsem_ByTargLocation{targpos}{replength}; zscore_sem_all];
        end
        
        % ==== put into output
        if length(RepeatFFmeanNORM_ByTargLocation{targpos})<replength;
            RepeatFFmeanNORM_ByTargLocation{targpos}{replength}=zscore_all_norm;
        else
            RepeatFFmeanNORM_ByTargLocation{targpos}{replength}=[RepeatFFmeanNORM_ByTargLocation{targpos}{replength}; zscore_all_norm];
        end
      
        % ==== put into output
        if length(RepeatFFsemNORM_ByTargLocation{targpos})<replength;
            RepeatFFsemNORM_ByTargLocation{targpos}{replength}=zscore_sem_all_norm;
        else
            RepeatFFsemNORM_ByTargLocation{targpos}{replength}=[RepeatFFsemNORM_ByTargLocation{targpos}{replength}; zscore_sem_all_norm];
        end
      
        
        
        
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


%% PLOT -- using non-normalized values of shift

lt_figure; hold on;

NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));
hsplots=[];

for i=1:NumTargLocations;
    
    hsplot=lt_subplot(NumTargLocations, 1, i);
    title(['target location: ' num2str(i)]);
    
    hold on;
    hsplots=[hsplots hsplot];
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
        % ====== PLOT INDIVIUDAL EXPERIMENTS
        x=1:ii;
        
        plot(x, RepeatFFmean_ByTargLocation{i}{ii}, '-ok');
        
        % ===== plot mean for this replength class
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
        Ymean=mean(RepeatFFmean_ByTargLocation{i}{ii});
        Ysem=lt_sem(RepeatFFmean_ByTargLocation{i}{ii});
        
        lt_plot(x, Ymean, {'LineStyle','-','Errors', Ysem, 'Color',[rand rand rand]});
        end
    end
    
    
    
end

    
% lt_plot_equalyaxis(hsplots);

lt_plot_equalaxes(hsplots,1,1);
% linkaxes(hsplots, 'x','y');

lt_subtitle('not-normalized');


%% PLOT - [non normalized] [separate by rep class] [OLD - DOESN'T WORK BECUASE SUBPLOTS ARE NOT WELL POSITIONED]


% lt_figure; hold on;
% NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation)); % number of rows
% 
% %    % number of columns
% % tmpfunc=@(x) cellfun(@isempty, x);
% % MaxNumRepLengths=cellfun(tmpfunc, RepeatFFmean_ByTargLocation, 'UniformOutput',0)
% % hsplots=[];
% 
% cc=1;
% for i=1:NumTargLocations;
%     
%     NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
%     
%     %     for ii=1:NumRepClasses;
%     
%     NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
%     
%     for ii=1:NumClassesOfRepeatLengths;
%         if isempty(RepeatFFmean_ByTargLocation{i}{ii})
%             continue
%         end
%         
%         hsplot=lt_subplot(NumTargLocations, NumRepClassesActual, cc);
%         title(['targ loc: ' num2str(i), '; rep class: ' num2str(ii)]);
%         
%         hold on;
%         hsplots=[hsplots hsplot];
%         
%         % ====== PLOT INDIVIUDAL EXPERIMENTS
%         x=1:ii;
%         
%         for k=1:size(RepeatFFmean_ByTargLocation{i}{ii},1);
%             lt_plot(x, RepeatFFmean_ByTargLocation{i}{ii}(k,:), {'Errors', RepeatFFsem_ByTargLocation{i}{ii}(k,:), 'LineStyle','-', 'MarkerFaceColor','none'});
%             if plot_text_names==1;
%                 string_to_plot=BirdExptTargNames{i}{ii}{k};
%                 
%                 text(x(end), RepeatFFmean_ByTargLocation{i}{ii}(k, end), string_to_plot);
%                 
%             end
%             
%         end
%         
%         
%         % ===== plot mean for this replength class
%         if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
%             
%             Ymean=mean(RepeatFFmean_ByTargLocation{i}{ii});
%             Ysem=lt_sem(RepeatFFmean_ByTargLocation{i}{ii});
%             
%             %             lt_plot(x, Ymean, {'LineStyle','-','Errors', Ysem, 'Color',[rand rand rand]});
%             lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','r'})
%         end
%         cc=cc+1;
%         
%         lt_plot_zeroline;
%         
%     end
%     %     end
% end
% 
% 
% % lt_plot_equalyaxis(hsplots);
% 
% lt_plot_equalaxes(hsplots,1,1);
% % linkaxes(hsplots, 'x','y');
% 
% 
% 
% % ++++++++++++++++++++++++++++++++ PLOT -- [normalized] [separating by class and targ loc]
% 
% 
% lt_figure; hold on;
% NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));
% 
% hsplots=[];
% 
% cc=1;
% for i=1:NumTargLocations;
%     
%     NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
%     
%     %     for ii=1:NumRepClasses;
%     
%     NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
%     
%     for ii=1:NumClassesOfRepeatLengths;
%         if isempty(RepeatFFmean_ByTargLocation{i}{ii})
%             continue
%         end
%         
%         hsplot=lt_subplot(NumTargLocations, NumRepClassesActual, cc);
%         title(['targ loc: ' num2str(i), '; rep class: ' num2str(ii)]);
%         
%         hold on;
%         hsplots=[hsplots hsplot];
%         
%         % ====== PLOT INDIVIUDAL EXPERIMENTS
%         x=1:ii;
%         for k=1:size(RepeatFFmean_ByTargLocation{i}{ii},1);
%             
%             % normalized data
%             zscore_Norm=RepeatFFmeanNORM_ByTargLocation{i}{ii}(k,:);
%             zscore_sem_Norm=RepeatFFsemNORM_ByTargLocation{i}{ii}(k,:);
%             
%             
%             % --- PLOT
%             lt_plot(x, zscore_Norm, {'Errors', zscore_sem_Norm, 'LineStyle','-', 'MarkerFaceColor','none'});
%             
%             if plot_text_names==1;
%                 string_to_plot=BirdExptTargNames{i}{ii}{k};
%                 
%                 text(x(end), RepeatFFmean_ByTargLocation{i}{ii}(k, end), string_to_plot);
%             end
%             
%         end
%         
%         
%         % ===== plot mean for this replength class
%         if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
%             
%             Ymean=nanmean(RepeatFFmeanNORM_ByTargLocation{i}{ii});
%             Ysem=lt_sem(RepeatFFsemNORM_ByTargLocation{i}{ii});
%             
%             %             lt_plot(x, Ymean, {'LineStyle','-','Errors', Ysem, 'Color',[rand rand rand]});
%             lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','r'})
%         end
%         cc=cc+1;
%       
%         lt_plot_zeroline;
%         
%     end
%     %     end
% end
% 
% 
% 
% lt_plot_equalaxes(hsplots,1,1);
% ylim([-0.4 1.4]);



%% PLOT - [non normalized] [separate by rep class]


count=1;
subplotrows=4;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation)); % number of rows

   % number of columns
hsplots=[];
cc=1;
for i=1:NumTargLocations;
    
    NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
    
    %     for ii=1:NumRepClasses;
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title(['targ loc: ' num2str(i), '; rep class: ' num2str(ii)]);
        
        hold on;
        hsplots=[hsplots hsplot];
        
        % ====== PLOT INDIVIUDAL EXPERIMENTS
        x=1:ii;
        
        for k=1:size(RepeatFFmean_ByTargLocation{i}{ii},1);
            lt_plot(x, RepeatFFmean_ByTargLocation{i}{ii}(k,:), {'Errors', RepeatFFsem_ByTargLocation{i}{ii}(k,:), 'LineStyle','-', 'MarkerFaceColor','none'});
            if plot_text_names==1;
                string_to_plot=BirdExptTargNames{i}{ii}{k};
                
                text(x(end), RepeatFFmean_ByTargLocation{i}{ii}(k, end), string_to_plot);
                
            end
            
        end
        
        
        % ===== plot mean for this replength class
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
            Ymean=mean(RepeatFFmean_ByTargLocation{i}{ii});
            Ysem=lt_sem(RepeatFFmean_ByTargLocation{i}{ii});
            
            %             lt_plot(x, Ymean, {'LineStyle','-','Errors', Ysem, 'Color',[rand rand rand]});
            lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','r'})
        end
        cc=cc+1;
        
        lt_plot_zeroline;
        
        % --- INDICATE TARGET LOCATION
        plot(i,0,'^r','MarkerSize',10);
        
    end
    %     end
end


% lt_plot_equalyaxis(hsplots);

lt_plot_equalaxes(hsplots,1,1);
% linkaxes(hsplots, 'x','y');




%% PLOT -- [normalized] [separating by class and targ loc]


count=1;
subplotrows=4;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));

hsplots=[];

cc=1;
for i=1:NumTargLocations;
    
    NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
    
    %     for ii=1:NumRepClasses;
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
[fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title(['targ loc: ' num2str(i), '; rep class: ' num2str(ii)]);
        
        hold on;
        hsplots=[hsplots hsplot];
        
        % ====== PLOT INDIVIUDAL EXPERIMENTS
        x=1:ii;
        for k=1:size(RepeatFFmean_ByTargLocation{i}{ii},1);
            
            % normalized data
            zscore_Norm=RepeatFFmeanNORM_ByTargLocation{i}{ii}(k,:);
            zscore_sem_Norm=RepeatFFsemNORM_ByTargLocation{i}{ii}(k,:);
            
            
            % --- PLOT
            lt_plot(x, zscore_Norm, {'Errors', zscore_sem_Norm, 'LineStyle','-', 'MarkerFaceColor','none'});
            
            if plot_text_names==1;
                string_to_plot=BirdExptTargNames{i}{ii}{k};
                
                text(x(end), RepeatFFmean_ByTargLocation{i}{ii}(k, end), string_to_plot);
            end
            
        end
        
        
        % ===== plot mean for this replength class
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
            Ymean=nanmean(RepeatFFmeanNORM_ByTargLocation{i}{ii});
            Ysem=lt_sem(RepeatFFsemNORM_ByTargLocation{i}{ii});
            
            %             lt_plot(x, Ymean, {'LineStyle','-','Errors', Ysem, 'Color',[rand rand rand]});
            lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','r'})
        end
        cc=cc+1;
      
        lt_plot_zeroline;
        
    end
    %     end
end



lt_plot_equalaxes(hsplots,1,1);
ylim([-0.4 1.4]);


%% ===== ONE PLOT FOR EACH TARG LOC, OVERLAY MEANS FOR EACH REP CLASS [normalized]
lt_figure; hold on;
NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));

hsplots=[];

for i=1:NumTargLocations;
    
    %     NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
%     hsplot=lt_subplot(NumTargLocations, 1, i);
    hsplot=lt_subplot(2,2, i);
    hold on;
    title(['targ loc: ' num2str(i)]);
    
    hsplots=[hsplots hsplot];
    
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
            Ymean=nanmean(RepeatFFmeanNORM_ByTargLocation{i}{ii});
            Ysem=lt_sem(RepeatFFsemNORM_ByTargLocation{i}{ii});
            
        else
            % if single, then just take single, no SEM
            Ymean=RepeatFFmeanNORM_ByTargLocation{i}{ii};
            Ysem=[];
        end
        
        if all(isnan(Ymean))
            continue
        end
        
        % ==== PLOT THIS REP CLASS
        x=1:ii;
        if isempty(Ysem);
            lt_plot(x, Ymean, {'LineStyle','-','MarkerFaceColor','none','Color',[rand rand rand]});
            
        else
            
            lt_plot(x, Ymean, {'LineStyle','-','Errors',Ysem, 'Color',[0.2 rand rand]});
        end
        
    end
        % ==== GLOBAL MEAN FOR THIS TARGET (might have different sample
        % sizes for at each syl location)
        
        Ymean_ALL=[];
        Ysem_ALL=[];
        for pos_compiling=1:length(RepeatFFmeanNORM_ByTargLocation{i});
            
            Y=[];
            for j=pos_compiling:length(RepeatFFmeanNORM_ByTargLocation{i});
                
                if isempty(RepeatFFmeanNORM_ByTargLocation{i}{j})
                    continue
                end
                
                y=RepeatFFmeanNORM_ByTargLocation{i}{j}(:,pos_compiling);
                Y=[Y; y];

            end
            
            if length(Y)<2
                break
            end
                
            Ymean_ALL(pos_compiling)=nanmean(Y);
            Ysem_ALL(pos_compiling)=lt_sem(Y);
        end
        
        lt_plot_zeroline;
        if isempty(Ymean_ALL) | all(isnan(Ymean_ALL))
            continue
        end
            
        
        lt_plot_bar(1:length(Ymean_ALL), Ymean_ALL, {'Errors',Ysem_ALL});
        line(xlim, [1 1], 'Color','k','LineStyle','--');
        
end



lt_plot_equalaxes(hsplots,1,1);
% ylim([-0.4 1.4]);

lt_subtitle('each expt normalized to targ [global mean across all rends]');



%% ===== ONE PLOT FOR EACH TARG LOC, OVERLAY MEANS FOR EACH REP CLASS [normalized]
lt_figure; hold on;
NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));

hsplots=[];

for i=1:NumTargLocations;
    
    %     NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
%     hsplot=lt_subplot(NumTargLocations, 1, i);
        hsplot=lt_subplot(2,2, i);
hold on;
    title(['targ loc: ' num2str(i)]);
    
    hsplots=[hsplots hsplot];
    
    Y_ALL={};
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
            Ymean=nanmean(RepeatFFmeanNORM_ByTargLocation{i}{ii});
            Ysem=lt_sem(RepeatFFsemNORM_ByTargLocation{i}{ii});
            
        else
            % if single, then just take single, no SEM
            Ymean=RepeatFFmeanNORM_ByTargLocation{i}{ii};
            Ysem=[];
        end
        
                if all(isnan(Ymean))
            continue
        end

        
        
        % ==== PLOT THIS REP CLASS
        x=1:ii;
        if isempty(Ysem);
            lt_plot(x, Ymean, {'LineStyle','-','MarkerFaceColor','none','Color',[rand rand rand]});
            
        else
            
            lt_plot(x, Ymean, {'LineStyle','-','Errors',Ysem, 'Color',[0.2 rand rand]});
        end
        
        Y_ALL{ii}=Ymean;
    end
    
        % ==== GLOBAL MEAN FOR THIS TARGET (might have different sample
        % sizes for at each syl location)
        
        Ymean_ALL=[];
        Ysem_ALL=[];
        for pos_compiling=1:length(Y_ALL);
            
            Y=[];
            for j=pos_compiling:length(Y_ALL);
                
                if isempty(Y_ALL{j})
                    continue
                end
                
                y=Y_ALL{j}(:,pos_compiling);
                Y=[Y; y];
               
            end
            
            if length(Y)<2
                break
            end
                
            Ymean_ALL(pos_compiling)=nanmean(Y);
            Ysem_ALL(pos_compiling)=lt_sem(Y);
        end
        
        lt_plot_zeroline;
                if isempty(Ymean_ALL) | all(isnan(Ymean_ALL))
            continue
        end

        
        lt_plot_bar(1:length(Ymean_ALL), Ymean_ALL, {'Errors',Ysem_ALL});
        line(xlim, [1 1], 'Color','k','LineStyle','--');
        
end



lt_plot_equalaxes(hsplots,1,1);
% ylim([-0.4 1.4]);

lt_subtitle('each expt normalized to targ, global mean across means');


    
%% ===== ONE PLOT FOR EACH TARG LOC, OVERLAY MEANS FOR EACH REP CLASS [not normalized] [golobal mean across rends]

lt_figure; hold on;
NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));

hsplots=[];

for i=1:NumTargLocations;
    
    %     NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
%     hsplot=lt_subplot(NumTargLocations, 1, i);
        hsplot=lt_subplot(2,2, i);
hold on;
    title(['targ loc: ' num2str(i)]);
    
    hsplots=[hsplots hsplot];
    
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
            Ymean=nanmean(RepeatFFmean_ByTargLocation{i}{ii});
            Ysem=lt_sem(RepeatFFsem_ByTargLocation{i}{ii});
            
        else
            % if single, then just take single, no SEM
            Ymean=RepeatFFmean_ByTargLocation{i}{ii};
            Ysem=[];
        end
        
        
        % ==== PLOT THIS REP CLASS
        x=1:ii;
        if isempty(Ysem);
            lt_plot(x, Ymean, {'LineStyle','-','MarkerFaceColor','none','Color',[rand rand rand]});
            
        else
            
            lt_plot(x, Ymean, {'LineStyle','-','Errors',Ysem, 'Color',[0.2 rand rand]});
        end
        
    end
        % ==== GLOBAL MEAN FOR THIS TARGET (might have different sample
        % sizes for at each syl location)
        
        Ymean_ALL=[];
        Ysem_ALL=[];
        for pos_compiling=1:length(RepeatFFmean_ByTargLocation{i});
            
            Y=[];
            for j=pos_compiling:length(RepeatFFmean_ByTargLocation{i});
                
                if isempty(RepeatFFmean_ByTargLocation{i}{j})
                    continue
                end
                
                y=RepeatFFmean_ByTargLocation{i}{j}(:,pos_compiling);
                Y=[Y; y];

            end
            
            if length(Y)<2
                break
            end
                
            Ymean_ALL(pos_compiling)=nanmean(Y);
            Ysem_ALL(pos_compiling)=lt_sem(Y);
        end
        
        lt_plot_zeroline;
                if isempty(Ymean_ALL) | all(isnan(Ymean_ALL))
            continue
        end

        
        lt_plot_bar(1:length(Ymean_ALL), Ymean_ALL, {'Errors',Ysem_ALL});
        line(xlim, [1 1], 'Color','k','LineStyle','--');
        
end



lt_plot_equalaxes(hsplots,1,1);
% ylim([-0.4 1.4]);

lt_subtitle('each expt normalized to targ [global mean across all rends]');


%% ===== ONE PLOT FOR EACH TARG LOC, OVERLAY MEANS FOR EACH REP CLASS [NOT normalized] [mean of means[
lt_figure; hold on;
NumTargLocations=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation));

hsplots=[];

for i=1:NumTargLocations;
    
    %     NumRepClassesActual=sum(~cellfun(@isempty, RepeatFFmean_ByTargLocation{i}));
    
    NumClassesOfRepeatLengths=length(RepeatFFmean_ByTargLocation{i});
    
%     hsplot=lt_subplot(NumTargLocations, 1, i);
    hsplot=lt_subplot(2,2, i);
    hold on;
    title(['targ loc: ' num2str(i)]);
    
    hsplots=[hsplots hsplot];
    
    Y_ALL={};
    for ii=1:NumClassesOfRepeatLengths;
        if isempty(RepeatFFmean_ByTargLocation{i}{ii})
            continue
        end
        
        if size(RepeatFFmean_ByTargLocation{i}{ii},1)>1;
            
            Ymean=nanmean(RepeatFFmean_ByTargLocation{i}{ii});
            Ysem=lt_sem(RepeatFFsem_ByTargLocation{i}{ii});
            
        else
            % if single, then just take single, no SEM
            Ymean=RepeatFFmean_ByTargLocation{i}{ii};
            Ysem=[];
        end
        
        
        % ==== PLOT THIS REP CLASS
        x=1:ii;
        if isempty(Ysem);
            lt_plot(x, Ymean, {'LineStyle','-','MarkerFaceColor','none','Color',[rand rand rand]});
            
        else
            
            lt_plot(x, Ymean, {'LineStyle','-','Errors',Ysem, 'Color',[0.2 rand rand]});
        end
        
        Y_ALL{ii}=Ymean;
    end
        % ==== GLOBAL MEAN FOR THIS TARGET (might have different sample
        % sizes for at each syl location)
        
        Ymean_ALL=[];
        Ysem_ALL=[];
        for pos_compiling=1:length(Y_ALL);
            
            Y=[];
            for j=pos_compiling:length(Y_ALL);
                
                if isempty(Y_ALL{j})
                    continue
                end
                
                y=Y_ALL{j}(:,pos_compiling);
                Y=[Y; y];
               
            end
            
            if length(Y)<2
                break
            end
                
            Ymean_ALL(pos_compiling)=nanmean(Y);
            Ysem_ALL(pos_compiling)=lt_sem(Y);
        end
        
        lt_plot_zeroline;
                if isempty(Ymean_ALL) | all(isnan(Ymean_ALL))
            continue
        end

        lt_plot_bar(1:length(Ymean_ALL), Ymean_ALL, {'Errors',Ysem_ALL});
        line(xlim, [1 1], 'Color','k','LineStyle','--');
        
end



lt_plot_equalaxes(hsplots,1,1);
% ylim([-0.4 1.4]);

lt_subtitle('no normalization, global mean across means');

    




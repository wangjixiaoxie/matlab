function SeqDepPitch_AcrossBirds = lt_seq_dep_pitch_ACROSSBIRDS_LearningRate(SeqDepPitch_AcrossBirds, ExtractRaw)
%% lt 5/26/17 - estimates learning rate by counting number of syls in a day, and getting hz/rend
% NOTE: CURRENTLY assumes that the target context/syl uses notenum 0 in
% evtaf online. 

% NOTE: only run this on single context training currently

lt_figure;
lt_plot_text(0, 0.5, 'NOTE: only do single context training expt... ', 'r');


NumBirds = length(SeqDepPitch_AcrossBirds.birds);


%% go through all expt and go through first 4 days of learning and extract renditions
if ExtractRaw==1
for i=1:NumBirds
   
    NumExpts = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:NumExpts
        
       
        % what are the learning days?
        numEmptyDays = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
        LearningDaysUsed = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;

        assert(all(LearningDaysUsed == ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.dayIndsUsed), 'asfasdf');
        assert(WNday1+numEmptyDays + 2 == LearningDaysUsed(1), 'asfsdaf');
        
        % ==== for each laerning day go to raw dat folder and extract num
        % renditions
        LearningDays = [WNday1+numEmptyDays:LearningDaysUsed(2)];
        assert(length(LearningDays) == 4, 'afasdf');
        
        FirstDay_date = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay; % first day of expt (including baseline)
        
        
        
        for j=1:length(LearningDays)
           dayInd = LearningDays(j);
           
           dayInd_date = datestr(datenum(FirstDay_date, 'ddmmmyyyy') + dayInd - 1, 'mmddyy');
           
           % -- go to folder with raw day
           cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir)
           cd ..
           cd ..
           
           % -- find day folder
           exptphrase = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.DayRawDat.phrase;
           
           tmp = dir([dayInd_date '_' exptphrase '_durWN*']);
           if length(tmp) ~=1
               tmp = dir([dayInd_date '_' exptphrase '*']);
           end
           assert(length(tmp) ==1, 'asfa');
           cd(tmp.name);
           
           %==============================================================
           % EXTRACT NUM RENDS THIS DAY
           if exist('CollectNumTrigsStruct.mat', 'file')
               % then skip, since already done for this day
               disp('SKIPPING DAY, alrady done');
               continue
           end
           
           
           CollectNumTrigsStruct = struct;
           [TotalNumTrigs, NotenumsUnique, NoteNumsAllVect] = lt_seq_dep_pitch_CollectNumTrigs;
           
%            if length(NotenumsUnique) >1

            % ================ SAVE
            CollectNumTrigsStruct.TotalNumTrigs = TotalNumTrigs;
            CollectNumTrigsStruct.NotenumsUnique = NotenumsUnique;
            CollectNumTrigsStruct.NoteNumsAllVect = NoteNumsAllVect;
            
            save('CollectNumTrigsStruct', 'CollectNumTrigsStruct');    
            
        end 
    end
end
end

%% ====== GO THRU ALL DAYS AND ESTIMATE LEARNING RATE

for i=1:NumBirds
   
    NumExpts = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:NumExpts
        
       
        % what are the learning days?
        numEmptyDays = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
        LearningDaysUsed = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.LearningMetric.dayIndsUsed;
        targsyl = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;

        assert(all(LearningDaysUsed == ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.dayIndsUsed), 'asfasdf');
        assert(WNday1+numEmptyDays + 2 == LearningDaysUsed(1), 'asfsdaf');
        
        % ==== for each laerning day go to raw dat folder and extract num
        % renditions
        LearningDays = [WNday1+numEmptyDays:LearningDaysUsed(2)];
        assert(length(LearningDays) == 4, 'afasdf');
        
        FirstDay_date = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay; % first day of expt (including baseline)
        
        
        
        for j=1:length(LearningDays)
           dayInd = LearningDays(j);
           
           dayInd_date = datestr(datenum(FirstDay_date, 'ddmmmyyyy') + dayInd - 1, 'mmddyy');
           
           % -- go to folder with raw day
           cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir)
           cd ..
           cd ..
           
           % -- find day folder
           exptphrase = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.DayRawDat.phrase;
           
           tmp = dir([dayInd_date '_' exptphrase '_durWN*']);
           if length(tmp) ~=1
               tmp = dir([dayInd_date '_' exptphrase '*']);
           end
           assert(length(tmp) ==1, 'asfa');
           cd(tmp.name);
           
           % =========== LOAD PREVIOUSLY SAVED DAT
           clear CollectNumTrigsStruct
           load('CollectNumTrigsStruct.mat'); 
           
           % ============ ESTIMATE NUM RENDS OF TARGET CONTEXT FOR THIS DAY
           numhits = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(targsyl).HIT_RATE_overdays.NumHits(dayInd);
           numtotal = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(targsyl).HIT_RATE_overdays.NumTotal(dayInd);
           
           % --- figure out false targ rate
           SylsUnique = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
           OffTarghits = [];
           for jj=1:length(SylsUnique)
               syl = SylsUnique{jj};
               hitstmp = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits(dayInd);
               
               OffTarghits = [OffTarghits hitstmp];
           end
           numhits_alltarg = sum(OffTarghits);
           
           
           NumHitsTotal = CollectNumTrigsStruct.TotalNumTrigs;
           
           NumRendsTotal = NumHitsTotal * (numhits/(numhits_alltarg)) ...
               * (numtotal/numhits); % scale by pr(on target) and pr is hit given is syl
% 
%            if isnan(NumRendsTotal)
%                keyboard
%            end
           
           disp([num2str(NumRendsTotal) ' rends; ' num2str(NumHitsTotal) ' hits']);
           
           % ================ OUTPUT
           SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.LEARNING_RATE_AT_TARG.NumRendsUnlabeledSongs(dayInd) = NumRendsTotal;
           
        end 
    end
end


%% ACCOUNT for offtarget hits by estimating hit rate using labeled songs.




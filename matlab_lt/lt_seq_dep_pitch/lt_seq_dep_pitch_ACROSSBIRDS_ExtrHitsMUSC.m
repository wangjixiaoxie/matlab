function lt_seq_dep_pitch_ACROSSBIRDS_ExtrHitsMUSC(SeqDepPitch_AcrossBirds)

%% LT 4/1/16 - extracts hit rate data for musc songs, because I was stupid and did not extract in early code
% extracts it from the raw dat struct

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% == go through all experiments, and if required, extract hit rate for musc songs
for i=1:NumBirds;
    birdname = SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % ===== skip if has musc data
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0;
            disp(['SKIPPED ' birdname ' - ' exptname ', not MUSC experiment']);
            continue
        end
        
        % ===== skip if already has hit rate for musc extracted
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl), 'HitRateStatus');
            
            disp(['SKIPPED ' birdname ' - ' exptname ', HitRateStatus already extracted!']);
            continue
        end
        
       
        %% ===== load params and alldaysplotlearning
        cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
               
        load Params
        load AllDays_PlotLearning
        load AllDays_RawDatStruct
        
        %% == COLLECT HIT RATE
        SylFields_Unique=Params.PlotLearning.SylFields_Unique;
        NumDays=length(AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).FFvals);
        
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            for k=1:NumDays;
                
                if isempty(AllDays_RawDatStruct{k})
                    continue
                end
                
                if ~isfield(AllDays_RawDatStruct{k}.data_MUSC, syl);
                    continue;
                end
                
                if isempty(AllDays_RawDatStruct{k}.data_MUSC.(syl));
                    continue
                end
                
                % collect data
                AllDays_PlotLearning.DataMatrix_MUSC.(syl).HitRateStatus{k}=...
                    cell2mat(AllDays_RawDatStruct{k}.data_MUSC.(syl)(:,9));
                
            end
        end
        
        %% ====== SAVE IN BIRD EXPT FOLDER
            save('AllDays_PlotLearning','AllDays_PlotLearning');
        disp(['DONE! - ' birdname '-' exptname]);
        
        lt_figure; hold on;
        lt_plot_annotation(1, 'HAVE TO RESTART lt_seq_dep_pitch', 'r');
        lt_plot_annotation(2, 'Reason: overwrote AllDaysPlotLearning, need to reload', 'k');
            
    end
end

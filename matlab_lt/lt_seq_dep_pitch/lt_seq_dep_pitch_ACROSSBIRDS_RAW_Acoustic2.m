function lt_seq_dep_pitch_ACROSSBIRDS_RAW_Acoustic2(ExperimentList)
% --- moves structure stats and DONE indicator to subfolder. does not touch
% MUSC struct stats.

currdir=pwd;
disp(' ');

for i=1:length(ExperimentList);
    
    if isempty(ExperimentList{i});
        continue
    end
    
    birdname=ExperimentList{i}{1};
    exptname=ExperimentList{i}{2};
    savedir=ExperimentList{i}{4};
    
    disp([ ' ---' birdname '-' exptname]);
    
    
    cd(savedir);
    
    
    if (0)
        % ==== backup old structs
        if exist('AllDays_StructStatsStruct.mat', 'file')==2
            % then file exists
            
            try cd BACKUP_StructStatsStruct_UsingPowerNotAmpl %^ test if folder exists
                cd ..
            catch err
                mkdir BACKUP_StructStatsStruct_UsingPowerNotAmpl
            end
            
            % -- move struct
            eval('!mv AllDays_StructStatsStruct.mat BACKUP_StructStatsStruct_UsingPowerNotAmpl/');
            disp('moved: AllDays_StructStatsStruct.mat');
            
            % -- move DONE indicators
            fname=dir('DONE_StructureStats_*');
            
            for j=1:length(fname)
                fn=fname(j).name;
                
                eval(['!mv ' fn ' BACKUP_StructStatsStruct_UsingPowerNotAmpl/']);
                disp(['moved: ' fn]);
            end
        end
    end
    
    
    
    % ======== restore backups
    % -- take from using amplitude
    if exist('BACKUP_StructStatsStruct', 'dir')==7
        % then folder exists
        
        % --- MOVE UP ONE FOLDER
        cd BACKUP_StructStatsStruct
        
        % -- 1) move struct
        eval('!mv AllDays_StructStatsStruct.mat ..');
        disp('moved: AllDays_StructStatsStruct.mat');
        
        % -- 2) move done indicators
        fname=dir('DONE_StructureStats_*');
        for j=1:length(fname)
            fn=fname(j).name;
            
            eval(['!mv ' fn ' ..']);
            disp(['moved: ' fn]);
        end
    end
end


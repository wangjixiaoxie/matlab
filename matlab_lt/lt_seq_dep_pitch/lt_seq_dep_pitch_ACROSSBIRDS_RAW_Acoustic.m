function lt_seq_dep_pitch_ACROSSBIRDS_RAW_Acoustic(ExperimentList)

%% STRATEGy
% 1) use regular expressions stuff to extract acoustic data [only works
% with stereotyped currently]

% 2) Get feature space using that data

% 3) Later on, using seq dep pitch, find position of each syl in the regexp string.

% 4) Save to birds own database so dont' have to redo


%%  

currdir=pwd;

for i=1:length(ExperimentList);
    
    if isempty(ExperimentList{i});
        continue
    end
    
    birdname=ExperimentList{i}{1};
    exptname=ExperimentList{i}{2};
    
    % --- 1) EXTRACT RAW DATA
    savedir=ExperimentList{i}{4};
    
    cd(savedir);
    load AllDays_RawDatStruct
    
    try 
        load AllDays_RegExpr
    catch err
        disp([birdname '-' exptname ' is missing AllDays_RegExpr - SKIPPING (perform PreProcess analysis first)']);
    end
    
    
    % ---- 2) FOR EACH MOTIF, AND THEN EACH SYL, EXTRACT BASELINE SOUND
    % DATA
    RegExpMotifList=ExperimentList{i}{7};
    
    for j=1:length(RegExpMotifList);
        Motif=RegExpMotifList{j};
        
        for k=1:length(Motif);
            syl=Motif{k};
        end
    end
    
            
    
    
    
    
end



function lt_extract_AllDaysPC(Params, date_range)

%% Run in bird directory

% e.g. input:

% Params.syl_list={'a','b'};
% Params.freq_range_list={[3000 4000], [1300 2200]};
% Params.pc_dur_list=[0.12, 0.14];
% Params.pc_harms_list=[1, 1];
% 
% Params.batch='batch.labeled.all';
% Params.experiment = 'CtxtDepPitch';
% 
% date_range=''; % e.g. {'20Apr2015','20May2015'}. leave blank ('') for all days



%%

NumSylTypes=length(Params.syl_list);
currdir=pwd;



%% ==== Check (and make if required) for save folder

foldername=[currdir '/extract_AllDaysPC_' Params.experiment];

try cd(foldername)
    cd(currdir);
catch err
    mkdir(foldername)
    disp('made new extract data folder !...');
end


%% ==== GO THROUGH EACH DAY AND EXTRACT PC

% 1) Collect dir information - find dirs that you want
MetadataStruct=lt_metadata_collect;

ListOfDirs=lt_metadata_find_dirs(MetadataStruct, Params.experiment, ...
    '', '', date_range, 1,2);


% 2) go through all dates and extract data

for i=1:length(ListOfDirs);
    
    cd(ListOfDirs(i).dirname);
    DatStruct=struct;
    
    for ii=1:NumSylTypes;
        syl=Params.syl_list{ii};
        freq_range=Params.freq_range_list{ii};
        pc_dur=Params.pc_dur_list(ii);
        pc_harms=Params.pc_harms_list(ii);
        
        % ======= extract stats and PC
        disp(' ');
        disp(['extracting stats and pitch contours, syl: ' syl])
        fvals = findwnote2tw_v4_LT(Params.batch, syl,'',-0.005,...
            freq_range, pc_dur*32000,1,1,'obs0',0);
        
        [PC_raw, pc_F, pc_T, ~]=jc_pitchcontourFV_LT(fvals,...
            1024,1020,1, freq_range(1),...
            freq_range(2),pc_harms,'obs0');
        disp('done ...');
        
        
        % ======= Compile useful stats into one structure
        numrends=length(fvals);
        
        for j=1:numrends;
            DatStruct.(syl)(j).filename=fvals(j).fn;% filename
            DatStruct.(syl)(j).datenum=fvals(j).datenum; % datenum
            DatStruct.(syl)(j).song_labels=fvals(j).lbl;    % labels
            DatStruct.(syl)(j).note_pos=fvals(j).NotePos; % note position
            DatStruct.(syl)(j).trig_status=fvals(j).TRIG;  % TRIG
            DatStruct.(syl)(j).catch_trial=fvals(j).CatchTrial;  % catch trial
            DatStruct.(syl)(j).catch_song=fvals(j).CatchSong;  % catch song
            DatStruct.(syl)(j).Pitch_contour=PC_raw(:,j); % pitch contour
            DatStruct.(syl)(j).pitch_contour_T=pc_T;
            DatStruct.(syl)(j).pitch_contour_F=pc_F;
        end
        
        
    end
    % ===== SAVE
    savedir=[foldername '/DatStruct_' ListOfDirs(i).date '_' ListOfDirs(i).condition '.mat'];
    save(savedir, 'DatStruct');
    
    params_savedir=[foldername '/Params_' ListOfDirs(i).date '_' ListOfDirs(i).condition '.mat'];
    save(params_savedir, 'Params');
    
    cd(currdir);
end


disp(' ');
disp('ALL DONE! ....');

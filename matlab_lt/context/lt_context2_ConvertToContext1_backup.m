function [AllSongsDataMatrix, Params_alldays]= lt_context2_ConvertToContext1(ALLDATSTRUCT, Params_global, syl)
%% LT 11/15/15 - convert output of lt_extract_CompilePC to format that is compatible with lt_context_PLOT
% inputs, e.g.:
% syl='b';

% TO DO: 
% 1) instead of all being notenum1 , use actual syl, thus allow alls yls to be in one dataset.

% OUTPUTS:
% note groups - 1 and 2, correspond to contexts that are in Params_global
% notenum = 1 (if running version with only one syl extracted)
%%

numdays = length(ALLDATSTRUCT);


%%
AllSongsDataMatrix=[];
filenames_cell={};

for i=1:numdays
    
    if isempty(ALLDATSTRUCT(i).data);
        continue
    end
    
    datmatrix_tmp=nan(length(ALLDATSTRUCT(i).data.(syl)), 9);
    
    ffvals=[ALLDATSTRUCT(i).data.(syl).FFvals];
    tvals=[ALLDATSTRUCT(i).data.(syl).Tvals_days];
    trigstat=[ALLDATSTRUCT(i).data.(syl).trig_status];
    
    filenames_cell=[filenames_cell {ALLDATSTRUCT(i).data.(syl).filename}];
    
    % condition - convert to note groups (1, 2, 3, ...)
    Condition_cell={ALLDATSTRUCT(i).data.(syl).condition};
    condition=nan(1,length(Condition_cell));
    for j=1:length(Params_global.CompilePC.all_conditions);
        cond=Params_global.CompilePC.all_conditions{j};
        
        
        condition(strcmp(Condition_cell, cond))=j;
    end
    
    datmatrix_tmp(:, 2)=tvals';
    datmatrix_tmp(:, 5)=ffvals';
    datmatrix_tmp(:, 7)=trigstat';
    datmatrix_tmp(:, 8)=condition';
    
    % ==== APPEND TO OUTPUT
    if isempty(AllSongsDataMatrix)
        AllSongsDataMatrix=datmatrix_tmp;
    else
        AllSongsDataMatrix=[AllSongsDataMatrix; datmatrix_tmp];
    end
end



%% -- Make legend
Params_global.ConvertToContext1.legend={'1_SongNum', '2_TimeOfSong_days','3_TriggerTime_IncludingPreBuffer','4_NoteNum','5_FF','6_Ampl','7_Hit','8_NoteGrougNum', '9_Experimental_Phase'};


%% make sure is in chron order
[~, inds]=sort(AllSongsDataMatrix(:, 2));
AllSongsDataMatrix=AllSongsDataMatrix(inds, :);



% unique songs
filenames_cell=filenames_cell(inds);
%     filenames_unique=unique(filenames_cell);
filenames_unique=unique(filenames_cell,'stable');

functmp=@(x) find(strcmp(filenames_unique, x));
songnum=cellfun(functmp, filenames_cell, 'UniformOutput', 0);
songnum=cell2mat(songnum);

AllSongsDataMatrix(:,1)=songnum';

%% Note number, set as 1 by default (since this is just one note)

AllSongsDataMatrix(:,4)=1;


%% Output params
Params_alldays=struct;

Params_alldays.firstday=Params_global.CompilePC.FirstDay;
Params_alldays.ConvertToContext1.all_conditions_in_order=Params_global.CompilePC.all_conditions;



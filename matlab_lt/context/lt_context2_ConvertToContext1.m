function [AllSongsDataMatrix, Params_alldays]= lt_context2_ConvertToContext1(ALLDATSTRUCT, Params_global)
%% LT 6/1/17 - If no notegroup data (i.e. no .ltrec2 file) then gives notegroup -1

%% LT 11/15/15 - convert output of lt_extract_CompilePC to format that is compatible with lt_context_PLOT
% inputs, e.g.:
% syl='b';

% TO DO: 
% 1) instead of all being notenum1 , use actual syl, thus allow alls yls to be in one dataset.

% OUTPUTS:
% note groups - 1 and 2, correspond to contexts that are in Params_global
% notenum = 1 (if running version with only one syl extracted)

% DEFINE NOTE GROUPS AND NOTE NUMBERS IN PARAMS


%%

numdays = length(ALLDATSTRUCT);

UseAutoNG = Params_global.ConvertToContext1.UseAutoNG;
if UseAutoNG==1
    Params_global.ConvertToContext1.NoteGroupNum_codes = {}; % delete entry, so no confusion later.
end

%%
AllSongsDataMatrix=[];
filenames_cell={};

% === perform once for each desired note:
num_notenums=length(Params_global.ConvertToContext1.NoteNum_codes)/2;

for n=1:num_notenums
    
    syl=Params_global.ConvertToContext1.NoteNum_codes{2*n-1};
    notenum=Params_global.ConvertToContext1.NoteNum_codes{2*n};
    
    for i=1:numdays
        
        if isempty(ALLDATSTRUCT(i).data);
            continue
        end
        
        % === isolate the dat struct, dependeing on if single syl or regexp
        if length(syl)==1;
            % then single syl
            DATSTRUCT_isolated=ALLDATSTRUCT(i).data.(syl);
        elseif length(syl)>1;
            % then is regexp
            DATSTRUCT_isolated=ALLDATSTRUCT(i).data_regexp.(syl);
        end
        
        
        datmatrix_output=nan(length(DATSTRUCT_isolated), 9); % for output
        
        ffvals=[DATSTRUCT_isolated.FFvals];
        tvals=[DATSTRUCT_isolated.Tvals_days];
        trigstat=[DATSTRUCT_isolated.trig_status];
        NoteNum_mat=notenum*ones(length(DATSTRUCT_isolated),1);
        
        filenames_cell=[filenames_cell {DATSTRUCT_isolated.filename}];
        
        
        % condition - convert to note groups (1, 2, 3, ...)
        if UseAutoNG==1
            condition_mat = [DATSTRUCT_isolated.NoteGroup];
        else
            Condition_cell={DATSTRUCT_isolated.condition};
            condition_mat=nan(1,length(Condition_cell));
            for j=1:length(Condition_cell)
                cond_current=Condition_cell{j};
                
                indtmp=find(strcmp(Params_global.ConvertToContext1.NoteGroupNum_codes, cond_current));
                
                if isempty(indtmp)
                    disp('PROBLEM - there are data with conditions not coded into Params_global.ConvertToContext1.NoteGroupNum_codes: GIVING nan for NGnum');
                    NGnum=nan;
                else
                    NGnum=Params_global.ConvertToContext1.NoteGroupNum_codes{indtmp+1};
                end
                
                condition_mat(j)=NGnum;
            end
        end
        
        if isempty(condition_mat)
            condition_mat = -1;
        end
        
        datmatrix_output(:, 2)=tvals';
        datmatrix_output(:, 5)=ffvals';
        datmatrix_output(:, 7)=trigstat';
        datmatrix_output(:, 8)=condition_mat';
        datmatrix_output(:, 4)=NoteNum_mat;
        
        % ==== APPEND TO OUTPUT
        if isempty(AllSongsDataMatrix)
            AllSongsDataMatrix=datmatrix_output;
        else
            AllSongsDataMatrix=[AllSongsDataMatrix; datmatrix_output];
        end
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


%% TROUBLESHOOTING - 
if (0) 
    AllSongsDataMatrix(:, [1 2 4 8])

rendnum=80;

ALLDATSTRUCT(1).data.b(rendnum).Tvals_days
ALLDATSTRUCT(1).data.b(rendnum).song_labels
ALLDATSTRUCT(1).data.b(rendnum).condition
end



%% Output params
Params_alldays=struct;

Params_alldays.firstday=Params_global.CompilePC.FirstDay;
Params_alldays.ConvertToContext1.all_conditions_in_order=Params_global.CompilePC.all_conditions;
Params_alldays.ConvertToContext1.NoteGroupNum_codes = Params_global.ConvertToContext1.NoteGroupNum_codes;
Params_alldays.ConvertToContext1.NoteNum_codes = Params_global.ConvertToContext1.NoteNum_codes;
Params_alldays.ConvertToContext1.UseAutoNG = Params_global.ConvertToContext1.UseAutoNG;


disp('DONE! ...');


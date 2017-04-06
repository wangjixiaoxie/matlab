function [ALLDATSTRUCT, Params_global]= lt_extract_CompilePC(plotON, Params_global, saveON)
% Go to extract folder made by lt_extract_AllDaysPC
%% PARAMS examples
% Params_global.CompilePC.PC_window_list={'b', [45 130], 'a', [365 460]}; % syl, value pairs
% Params_global.CompilePC.FirstDay='';
% Params_global.CompilePC.LastDay='';

% Regular expressions - first calculates FF etc, then performs regular
% expressions
% Params_global.CompilePC.regexp_list={'dcc(b)', 'bcc(b)'}; % e.g. {'dcc(b)', 'ab+(g)'} : dcc(b) means match dccb, and give me ind of b in parantheses.  ab+g means match ab... (anly length repeat), then g. give me ind of g
% Params_global.CompilePC.regexp_fieldnames={'dccB','bccB'}; % whatever
% want to call structure field (if this cell array not defined, then will
% attempt to use the regexp names.


%% go to folder (e.g.extract_AllDaysPC_CtxtDepPitch) where DatStruct and Params are saved - this compiles (1 - extract seq-dep-syl, and 2, calc FF)
TargStr='DatStruct*';
Dates_first='';
Dates_last='';
FileType='mat';
make_batch=0;

% get cell of dat structs in chrono order
[~, cell_of_structnames, Dates_first, Dates_last, arrayofdateinds]=lt_write_all_folder_contents_to_batch_v2(...
    TargStr, Params_global.CompilePC.FirstDay, Params_global.CompilePC.LastDay, FileType, make_batch);

Params_global.CompilePC.FirstDay=Dates_first;
Params_global.CompilePC.LastDay=Dates_last;


%% compile all data into one structure
ALLDATSTRUCT=struct;
unique_days=unique(arrayofdateinds);

for day=unique_days;
    inds=find(arrayofdateinds==day);
    
    ALLDATSTRUCT(day).data=[];
    % extract all dat structs for this day:
    for j=inds;
        load(cell_of_structnames{j});
        
        % === what was condition for this dat struct?
        uscores=strfind(cell_of_structnames{j},'_');
        if length(uscores)>1;
            % then condition is defined
            suffix=strfind(cell_of_structnames{j}, '.mat');
            condition=cell_of_structnames{j}(uscores(2)+1:suffix-1);
        else
            % then no condition defined
            condition = '';
        end
        
        syllist=fieldnames(DatStruct);
        % == for this condition, add condition information
        for k=1:length(syllist);
            syl=syllist{k};
            
            [DatStruct.(syl).condition]=deal(condition);
        end
        
        % === add to this day's structure
        for k=1:length(syllist);
            syl=syllist{k};
            
            if ~isfield(ALLDATSTRUCT(day).data, syl);
                ALLDATSTRUCT(day).data.(syl)=DatStruct.(syl);
            else
                ALLDATSTRUCT(day).data.(syl)=[ALLDATSTRUCT(day).data.(syl) DatStruct.(syl)];
            end
        end
        
        % ==================== PARAMS for this day
        params_fn=['Params' cell_of_structnames{j}(10:end)];
        load(params_fn);
        
        Params_global.ind_days(day)=Params;
        
    end
    
end

%% get some params

numdays=length(ALLDATSTRUCT);

% -- list of all syls used
allsyls={};
for i=1:numdays;
    try
        allsyls=[allsyls Params_global.ind_days(i).syl_list];
    catch err
    end
end
allsyls=unique(allsyls);
disp(['Found syllables (across all days): ' allsyls]);

Params_global.CompilePC.all_single_syls=allsyls;


% --- List of all conditions that exist
all_conditions={};
for i=1:numdays;
    try
        all_conditions=[all_conditions ALLDATSTRUCT(i).data.(allsyls{1}).condition];
    catch err
    end
end
all_conditions=unique(all_conditions);
disp(['Found conditions (across all days): ' all_conditions]);

Params_global.CompilePC.all_conditions=all_conditions;




%% Calculate FF from pitch contour
numdays=length(ALLDATSTRUCT);

% ===================== Plot all pitch contours
if plotON==1;
    
    for k=1:length(allsyls);
        syl=allsyls{k};
        
        count=1;
        SubplotsPerFig=10;
        subplotrows=5;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        
        
        
        hplotlist=[];
        for i=1:numdays
            
            if isempty(ALLDATSTRUCT(i).data);
                continue;
            end
            
            PC_mat=[ALLDATSTRUCT(i).data.(syl).Pitch_contour];
            
            [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
            plot(PC_mat);
            
            hplotlist=[hplotlist hsplot];
            
            axis tight
            title(['day' num2str(i)]);
            
            % plot time windows if they exist
            try
                ind=find(strcmp(Params_global.CompilePC.PC_window_list, syl))+1;
                timewindow=Params_global.CompilePC.PC_window_list{ind};
                
                line([timewindow(1) timewindow(1)], ylim, 'LineWidth',2,'Color','k');
                line([timewindow(2) timewindow(2)], ylim, 'LineWidth',2,'Color','k');
                
            catch err
            end
            
        end
        
        linkaxes(hplotlist, 'xy');
        lt_subtitle(syl);
        
    end
end



% ====================== Extract FF information
for i=1:numdays
    if isempty(ALLDATSTRUCT(i).data);
        continue,
    end
    
    for ii=1:length(allsyls);
        syl=allsyls{ii};
        
        
        PC_mat=[ALLDATSTRUCT(i).data.(syl).Pitch_contour];
        
        ind=find(strcmp(Params_global.CompilePC.PC_window_list, syl))+1;
        timewindow=Params_global.CompilePC.PC_window_list{ind};
        
        PC_mat_sliced=PC_mat(timewindow(1):timewindow(2),:);
        FFvals=mean(PC_mat_sliced,1);
        
        c=num2cell(FFvals);
        [ALLDATSTRUCT(i).data.(syl).FFvals]=deal(c{:});
        
    end
end

disp('FFvals extracted');


%% convert times to times in days rel to day 1
% ====================== Extract FF information
for i=1:numdays
    if isempty(ALLDATSTRUCT(i).data);
        continue,
    end
    
    for ii=1:length(allsyls);
        syl=allsyls{ii};
        
        
        datenum_list=[ALLDATSTRUCT(i).data.(syl).datenum];
        
        eventtimes=lt_convert_EventTimes_to_RelTimes(Params_global.CompilePC.FirstDay, datenum_list);
        c=num2cell(eventtimes.FinalValue);
        [ALLDATSTRUCT(i).data.(syl).Tvals_days]=deal(c{:});
        
    end
end


%% Sequence Filter if desired

if isfield(Params_global.CompilePC, 'regexp_list');
    % ===== for each regexpt - go through and find all cases where it exists in
    % song labels.  if the case has data, record down the ind and ff of that
    % data.
    Params_global.CompilePC.regexp_list_singlesyls={};
    
    for r=1:length(Params_global.CompilePC.regexp_list);
        regstring=Params_global.CompilePC.regexp_list{r};
        
        tmp1=findstr(regstring, '(');
        tmp2=findstr(regstring, ')');
        
        if tmp2-tmp1~=2;
            disp('PROBLEM - this code only works if the token (in parantheses) is a single syl');
        end
        
        if isempty(tmp1);
            disp('PROBLEM - a regexp does not have parantheses - need token to tell me which syl to calc pitch');
        end
        Singlesyl=regstring(tmp1+1:tmp2-1);
        Params_global.CompilePC.regexp_list_singlesyls{r}=Singlesyl;
        
        
        % +==================== EXTRACT DATA
        for i=1:numdays
            
            if isempty(ALLDATSTRUCT(i).data);
                ALLDATSTRUCT(i).data_regexp=[];
                continue;
            end
            
            numrends=length(ALLDATSTRUCT(i).data.(Singlesyl));
            
            % start field
            if isfield(Params_global.CompilePC, 'regexp_fieldnames');
                regstring_field=Params_global.CompilePC.regexp_fieldnames{r};
            else
                regstring_field= [regstring(1:tmp1-1) '_' regstring(tmp1+1) '_' regstring(tmp2+1:end)];  % replace parantheses with underscores
            end
            ALLDATSTRUCT(i).data_regexp.(regstring_field)=[];
            
            for ii=1:numrends;
                
                % for each rendition, ask whether it is in this sequence context
                sng_label=ALLDATSTRUCT(i).data.(Singlesyl)(ii).song_labels;
                note_pos=ALLDATSTRUCT(i).data.(Singlesyl)(ii).note_pos;
                
                
                [StartInd, MatchStr, tokenExtents]=regexp(sng_label, regstring, 'start', 'match', 'tokenExtents');
                
                % Is this rendition in this context?
                tmp=cell2mat(tokenExtents);
                tokenInds=tmp(1:2:end);
                
                if any(tokenInds==note_pos);
                    
                    % ----- REND IS IN CONTEXT - save it
                    TMPSTRUCT=struct;
                    
                    TMPSTRUCT.FFvals=ALLDATSTRUCT(i).data.(Singlesyl)(ii).FFvals;
                    TMPSTRUCT.Tvals_days=ALLDATSTRUCT(i).data.(Singlesyl)(ii).Tvals_days;
                    TMPSTRUCT.filename=ALLDATSTRUCT(i).data.(Singlesyl)(ii).filename;
                    TMPSTRUCT.trig_status=ALLDATSTRUCT(i).data.(Singlesyl)(ii).trig_status;
                    TMPSTRUCT.condition=ALLDATSTRUCT(i).data.(Singlesyl)(ii).condition;
                    TMPSTRUCT.NoteGroup = ALLDATSTRUCT(i).data.(Singlesyl)(ii).NoteGroup;
                    TMPSTRUCT.NotesInNoteGroup = ALLDATSTRUCT(i).data.(Singlesyl)(ii).NotesInNoteGroup;
                    
                    TMPSTRUCT.match_string=MatchStr{tokenInds==note_pos};
                    TMPSTRUCT.match_string_startind=StartInd(tokenInds==note_pos);
                    TMPSTRUCT.IND_SingleSylData=ii;
                    TMPSTRUCT.regexp_string=regstring;
                    TMPSTRUCT.SingleSyl=Singlesyl;
                    
                    
                    if isempty(ALLDATSTRUCT(i).data_regexp.(regstring_field));
                        ALLDATSTRUCT(i).data_regexp.(regstring_field)=TMPSTRUCT;
                    else
                        ALLDATSTRUCT(i).data_regexp.(regstring_field) = [ALLDATSTRUCT(i).data_regexp.(regstring_field) TMPSTRUCT];
                    end
                    
                end
            end
        end
        
        
        
    end
    
end
%% TROUBLESHOOTING regexp
% match reg exp filtered with original - make sure all same
if (0)
    % hand enter:
    day=4;
    regexp_field='dcc_b_';
    
    notes_before_targ=3;
    notes_post_targ=0;
    
    % run:
    numrends=length(ALLDATSTRUCT(day).data_regexp.(regexp_field));
    
    for k=1:numrends
        
        ffval_2=ALLDATSTRUCT(day).data_regexp.(regexp_field)(k).ffval;
        single_syl=ALLDATSTRUCT(day).data_regexp.(regexp_field)(k).SingleSyl;
        ind_orig=ALLDATSTRUCT(day).data_regexp.(regexp_field)(k).IND_SingleSylData;
        match_str_pos=ALLDATSTRUCT(day).data_regexp.(regexp_field)(k).match_string_startind;
        match_str=ALLDATSTRUCT(day).data_regexp.(regexp_field)(k).match_string;
        
        
        ffval_1=ALLDATSTRUCT(day).data.(single_syl)(ind_orig).FFvals;
        tmp1=ALLDATSTRUCT(day).data.(single_syl)(ind_orig).note_pos-notes_before_targ;
        tmp2=ALLDATSTRUCT(day).data.(single_syl)(ind_orig).note_pos+notes_post_targ;
        
        motif_orig=ALLDATSTRUCT(day).data.(single_syl)(ind_orig).song_labels(tmp1:tmp2);
        
        % display comparision of original and new to user
        
        disp([num2str(ffval_1) ' --- ' num2str(ffval_2)]);
        disp([match_str ' --- ' motif_orig]);
        disp(' ');
        
    end
end




%% SAVE
if saveON==1;
    try cd COMPILED
        
    catch err
        mkdir COMPILED
        cd COMPILED
    end
    
    save('ALLDATSTRUCT','ALLDATSTRUCT', '-v7.3');
    save('Params_global','Params_global');
else
    disp('DID NOT SAVE!....');
end







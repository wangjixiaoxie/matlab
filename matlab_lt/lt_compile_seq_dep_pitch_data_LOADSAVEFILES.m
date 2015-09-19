function SaveName=lt_compile_seq_dep_pitch_data_LOADSAVEFILES(input1,FirstDay,LastDay);
% RUN IN FOLDER CONTAINING THE saved files saved by
% lt_compile_seq_dep_pitch_data or lt_compile_seq_dep_pitch_data_SEQFILTER.
% This program opens all files and puts in order into one structure (based
% on date).

% input1=1; % will grab all files that exist and start structure at
% earliest day

% or:

% input1=0; % will have to specifify first and last day
% FirstDay='11Oct2014';
% LastDay='02Nov2014';

%%
% COLLECT ALL DATA into batch
batch=lt_write_all_folder_contents_to_batch(0);


% GO LINE BY LINE AND GET DATE AND FILE NAME
fid=fopen(batch);
c=1;
nextline=fgetl(fid);

while nextline~=-1; % i.e while files still exist
    filename{c}=nextline;
    filedate{c}=nextline(12:20);
    
    nextline=fgetl(fid);
    c=c+1;
end


% WHAT IS INDEX OF FIRST AND LAST DAY?
if input1==0; % then will use specified dates
    
    FirstDateNum=datenum(FirstDay,'ddmmmyyyy');
    LastDateNum=datenum(LastDay,'ddmmmyyyy');
elseif input1==1; % will figure out all dates
    
    FirstDateNum=min(datenum(filedate,'ddmmmyyyy'));
    LastDateNum=max(datenum(filedate,'ddmmmyyyy'));
    
    FirstDay=datestr(FirstDateNum,'ddmmmyyyy');
    LastDay=datestr(LastDateNum,'ddmmmyyyy');
end

% OPEN EACH FILE IN AND PLACE IN CORRECT POSITION (by date) IN LARGER
% STRUCTURE
AllDays_compiled_seqdep_pitch=cell(LastDateNum-FirstDateNum+1,1); % one cell for every desired day.
for i=1:length(filename);
    X=load(filename{i});
    
    ii=datenum(filedate{i},'ddmmmyyyy')-FirstDateNum+1; % date position
    
    if ii>0 && ii<=LastDateNum-FirstDateNum+1; % only compile data if it is within range of desired dates (only matters if I specified input dates)
    AllDays_compiled_seqdep_pitch{ii}=X.compiled_seqdep_pitch; % slot data into correct position
    end
    
end


%% SAVE

SaveName=['AllDays_Compiled_' FirstDay '_to_' LastDay];

try
    cd('AllDays_Compiled');
catch err
    mkdir('AllDays_Compiled');
    cd('AllDays_Compiled');
end

save(SaveName,'AllDays_compiled_seqdep_pitch','-v7.3');

cd ..

disp('DONE!');


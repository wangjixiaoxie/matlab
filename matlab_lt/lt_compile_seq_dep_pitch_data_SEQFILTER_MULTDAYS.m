function lt_compile_seq_dep_pitch_data_SEQFILTER_MULTDAYS(all_days,FirstDay,LastDay,SeqPreList,SylTargList,SeqPostList)
%% LT 11/18/14 - Load data structures containing single syls (e.g. a, b, c), and filters out sequence specific stuff (e.g. a[b])
% Run in the folder containing the data structures (e.g. % Run in folder: 
% e.g.
% /bluejay3/lucas/birds/pu11wh87/compile_seq_dep_pitch_data_SeqDepPitchShift,
% containing structures:
% e.g. DataStruct_02Nov2014);

%% PARAMETERS
% all_days = 1 (Run on all data in folder); =0 (Use the dates specified
% below). If 1, then doesn't matter what I enter for days argumemtns.
% FirstDay='11Oct2014'; %to filter out days I don't want
% LastDay='02Nov2014';

% WHAT SEQUENCES DO I WANT?
% SeqPreList={'ac','ab'}; % format: 1st elem of all three lists should
% combine. in this case the two sequences are ac[c]b and ab[b].
% SylTargList={'c','b'}; % these must already have raw data compiled above
% SeqPostList={'b',''};



%% LOAD SINGLE SYL DATA STRUCTURES

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
if all_days==0; % then will use specified dates
    
    FirstDateNum=datenum(FirstDay,'ddmmmyyyy');
    LastDateNum=datenum(LastDay,'ddmmmyyyy');
elseif all_days==1; % will figure out all dates
    
    FirstDateNum=min(datenum(filedate,'ddmmmyyyy'));
    LastDateNum=max(datenum(filedate,'ddmmmyyyy'));
    
    FirstDay=datestr(FirstDateNum,'ddmmmyyyy');
    LastDay=datestr(LastDateNum,'ddmmmyyyy');
end

% OPEN EACH FILE IN AND PLACE IN CORRECT POSITION (by date) IN LARGER
% STRUCTURE
for i=1:length(filename);
    ii=datenum(filedate{i},'ddmmmyyyy')-FirstDateNum+1; % date position (to check if I want this structure)

    if ii>0 && ii<=LastDateNum-FirstDateNum+1; % only compile data if it is within range of desired dates (only matters if I specified input dates)
        X=load(filename{i});
        lt_compile_seq_dep_pitch_data_SEQFILTER(X.compiled_seqdep_pitch,SeqPreList,SylTargList,SeqPostList); % Run function on the data structure just loaded - will save to a specific folder: savedir=['/bluejay' bluejaynum '/lucas/birds/' birdname '/compile_seq_dep_pitch_data_' phrase '/SEQFILTER'];
    end
end

disp('DONE!');


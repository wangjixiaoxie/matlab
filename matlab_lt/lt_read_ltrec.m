function RecDatStruct=lt_read_ltrec(fname_ltrec)
%% LT 3/7/15 - Load .ltrec file and output data as structure
% .ltrec file tracks all triggers. keeps track of ff, hit status,
% threshold, etc.  Useful for version of evtaf with autothresholder,
% modified by LT.

% Instructions: run in same directory containing the .ltrec file.  This will
% save the output structure in same directory

% inputs example:
% fname_ltrec='rec_file.ltrec';


%% Extract rec data into cell array
disp('Extracting .ltrec data into data structure...');

fid1=fopen(fname_ltrec); 
RecDatCell=textscan(fid1,'%f %s %s %s %s %s %s');


%% Put different rec data into separate fields of output struct
RecDatStruct=struct;

% FF estimates (hits (catch and notcatch) and escapes
RecDatStruct.FF_AllDetects=RecDatCell{1};

% Date/Time of trigger
TrigTimes_mat=cell2mat(RecDatCell{2});
RecDatStruct.Datenum_msresolution=datenum(TrigTimes_mat,'ddmmyy_HHMMSS.FFF'); % resolution of ms.

% Hit status
RecDatStruct.HitStatus=RecDatCell{3}; % either 1) hit trig, 2) escape, 3) hit catch trial, or 4) hit catch song

% Filename of song (.cbin) during trigger
for i=1:length(RecDatCell{4});
    tmp=strfind(RecDatCell{4}{i},'\'); % want just cbin name, without path
    RecDatStruct.SongFname_cbin{i}=RecDatCell{4}{i}(tmp(end)+1:end);
end

% Time of trigger within song cbin file (in ms)
for i=1:length(RecDatCell{5});
    trigtime=RecDatCell{5}{i}(9:end-2); % in ms
    trigtime=str2num(trigtime);
    
    RecDatStruct.TrigTimes_WithinCbin(i)=trigtime;
end
    
% BELOW: USEFUL IF USED AUTO THRESHOLD UPDATER -----------------------------

% FF threshold applied to each trial 
for i=1:length(RecDatCell{6});
    ffthresh=RecDatCell{6}{i}(9:end);
    ffthresh=str2num(ffthresh);
    
    RecDatStruct.FFThresh(i)=ffthresh;
end

% Desired hit rate and how many previous renditions to avearge over to
% determine hit rate.
for i=1:length(RecDatCell{7});
    
    ind=strfind(RecDatCell{7}{i},'NumRends');
    
    desiredhr=RecDatCell{7}{i}(10:ind-1); % desired hit rate, in percent
    desiredhr=str2num(desiredhr);
    desiredhr=desiredhr/100; % convert back to fraction
    
    numrends=RecDatCell{7}{i}(ind+8:end); % number of renditions
    numrends=str2num(numrends);
    
    RecDatStruct.DesiredHitRate(i)=desiredhr;
    RecDatStruct.NumRendsToCalcHitRate(i)=numrends;
end


%% SAVE

savename=[fname_ltrec(1:end-6) '_structure']; % take name string excluding '.ltrec' part
save(savename, 'RecDatStruct');

disp('Done!')



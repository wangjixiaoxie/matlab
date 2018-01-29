function TrigInfo=lt_batchsong_triginfo(batch)
%% LT 9/14/15 - quickly goes thru all songs and extracts trig info. no need for offline templates, etc
% CURRENTLY only works with evtafv4 (date from filename, and trig format)

% OUTPUTS:
% NumFB_ActualTrigs % is actual trigs (i.e. not catch trial (and I beleive
% not song), and must pass thresholds

%% RUN

% ================ INITIATE OUTPUTS
SongDatenum_all=[];
NumFB_ActualTrigs_all=[];
SongFilename_all={};

% ================= GO THRU BATCH
fid1=fopen(batch,'r');
rec_fn=fgetl(fid1);


while ischar(rec_fn)
    
    % --------------------- LOAD
rec_fn=[rec_fn(1:end-5) '.rec'];
% rec_struct=readrecf(rec_fn);
rec_struct=readrecf_LT_evtafv4(rec_fn); 


% ---------------- FIND NUM TRIGS
num_FB_actualtrigs=sum(~cellfun(@isempty,strfind(rec_struct.pbname, 'FB')));
num_FB_catch=sum(~cellfun(@isempty,strfind(rec_struct.pbname, 'catch')));
trigcatch = ~cellfun(@isempty,strfind(rec_struct.pbname, 'catch'));

% --------------- TIMES OF TRIGGERS
trigtimes = rec_struct.ttimes;

% ---------------- FILE DATENUM
uscores=strfind(rec_fn, '_');
song_datenum=rec_fn(uscores(end-1)+1:uscores(end)+6);
song_datenum=datenum(song_datenum, 'ddmmyy_HHMMSS');

% ==== collect all
SongDatenum_all=[SongDatenum_all song_datenum];
NumFB_ActualTrigs_all=[NumFB_ActualTrigs_all num_FB_actualtrigs];
SongFilename_all=[SongFilename_all [rec_fn(1:end-4) '.cbin']];
    
rec_fn=fgetl(fid1);
end

% ===================== OUTPUT
TrigInfo.SongDatenum=SongDatenum_all;
TrigInfo.NumFB_ActualTrigs=NumFB_ActualTrigs_all;
TrigInfo.SongFilename=SongFilename_all;
TrigInfo.trigtimes = trigtimes;
TrigInfo.trigcatch = trigcatch;

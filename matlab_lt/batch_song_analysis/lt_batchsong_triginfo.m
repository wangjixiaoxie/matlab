function TrigInfo=lt_batchsong_triginfo(batch)
%% LT 9/14/15 - quickly goes thru all songs and extracts trig info. no need for offline templates, etc
% Note: not checked thoroughly that online trig info is accurate.
% only works with evtafv4 (date from filename, and trig format)

% OUTPUTS:
% NumFB_ActualTrigs % is actual trigs (i.e. not catch trial (and I beleive
% not song), and must pass thresholds

fid1=fopen(batch,'r');
rec_fn=fgetl(fid1);

SongDatenum_all=[];
NumFB_ActualTrigs_all=[];
SongFilename_all={};

while isstr(rec_fn);
    
rec_fn=[rec_fn(1:end-5) '.rec'];
rec_struct=readrecf(rec_fn);

num_FB_actualtrigs=sum(~cellfun(@isempty,strfind(rec_struct.pbname, 'FB')));

uscores=strfind(rec_fn, '_');
song_datenum=rec_fn(uscores(end-1)+1:uscores(end)+6);
song_datenum=datenum(song_datenum, 'ddmmyy_HHMMSS');

% ==== collect all
SongDatenum_all=[SongDatenum_all song_datenum];
NumFB_ActualTrigs_all=[NumFB_ActualTrigs_all num_FB_actualtrigs];
SongFilename_all=[SongFilename_all [rec_fn(1:end-4) '.cbin']];
    
rec_fn=fgetl(fid1);
end

TrigInfo.SongDatenum=SongDatenum_all;
TrigInfo.NumFB_ActualTrigs=NumFB_ActualTrigs_all;
TrigInfo.SongFilename=SongFilename_all;

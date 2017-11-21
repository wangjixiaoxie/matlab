%%  TO MAKE BOS FILES (with digitual pulse for syl onset/offset

songfname = 'pu69wh78_031117_102806.9552.cbin'; % shoudl have notmat also
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

% ######################## convert to wavfile
% -- chan 1 (song)
[dat, fs] = ReadCbinFile(songfname);

dat = dat/max(dat);

% --- bandpass filter
filter_type='hanningfirff';
F_low  = 500;
F_high = 8000;
dat=bandpass(dat,fs,F_low,F_high,filter_type);
datbeforesmth = dat;


% -- chan 2 (onsets and offsets)
notdat = load([songfname '.not.mat']);

datdig = zeros(length(dat),1);

for i=1:length(notdat.onsets)
    ons = notdat.onsets(i);
    off = notdat.offsets(i);
    
    % convert from ms to samps
    
    ons = round(fs*(ons/1000));
    off = round(fs*(off/1000));
    
    datdig(ons:off) = 1;
end

% ######################## exponential falloff from sound onset and offset
rolltime = 0.2;
expsize = round(rolltime*fs);

% -- roll on
smth = exp(-((expsize:-1:1)-1)/(expsize/5));
dat(1:expsize) = dat(1:expsize).*smth';

% --- roll off
smth = exp(-((1:expsize)-1)/(expsize/5));
dat(end-expsize+1:end) = dat(end-expsize+1:end).*smth';


figure; hold on;
plot(datbeforesmth, 'k');
plot(dat, 'c');
plot(datdig, 'r');


% ############# combine song and dig signal
fnameout = [songfname '_DigOnsOff.wav'];
fnameout = [OutDir '/' fnameout];

wavwrite([dat datdig], fs, 16, fnameout);

% ############# MAKE REVERSE VERSION
dat_rev = flipud(dat);
datdig_rev = flipud(datdig);

fnameout = [songfname '_DigOnsOff_REV.wav'];
fnameout = [OutDir '/' fnameout];
wavwrite([dat_rev datdig_rev], fs, 16, fnameout);

% ########### plot
figure; hold on;
subplot(211); hold on;
title('forward');
plot(datbeforesmth, 'k');
plot(dat, 'c');
plot(datdig, 'r');

subplot(212); hold on;
title('rev');
plot(dat_rev, 'c');
plot(datdig_rev, 'r');


%%   script for getting song names and syl positions, given only part of song name and the motif number
% songlist = {'9.124', '46.203', '06.290'}; % partially enter strings for song names
% motifnums = [4 2 4]; % rendition of the motif within the song
% motifname = 'jbh'; % motif
% sylnuminmotif = 3; % which syl do you want in the motif?
songlist = {'46.203'}; % partially enter strings for song names
motifnums = [4]; % rendition of the motif within the song
motifname = 'abh'; % motif
sylnuminmotif = 3; % which syl do you want in the motif?

disp('--')

songsout = {};
sylposout = [];
for i=1:length(songlist)
    tmp = dir(['*' songlist{i} '*.cbin']);
    disp(tmp(1).name);
    assert(length(tmp)==1, 'asdfasd');
    
    % -- find location of desired syl
    load([tmp(1).name '.not.mat']);
    
    inds = strfind(labels, motifname);
    
    sylpos = inds(motifnums(i))+sylnuminmotif-1;
    
    disp(['POSITION: ' num2str(sylpos)])
    disp(labels(sylpos))
    disp(labels(sylpos-2:sylpos+2))
    
    
    songsout = [songsout tmp(1).name];
    sylposout = [sylposout sylpos];
end

songsout
sylposout
%% ================== splicing syllables into a BOS file
clear all; close all;
songfname = '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin'; % --- load cbin file, splice syllables, then save as wav file
inds_abh = [18 34 49]; % inds within the song (not.mat)
inds_jbh = [28 43 58];
sylstoremove = [inds_abh inds_jbh];

% --- abh
replfiles_abh_hi = {'pu69wh78_101117_223206.290.cbin', ...
    'pu69wh78_101117_222646.203.cbin', 'pu69wh78_101117_203411.123.cbin'};
notes_abh_hi = [38 48 17];

replfiles_abh_lo = { 'pu69wh78_101117_114745.3.cbin', ...
    'pu69wh78_101117_114457.188.cbin', 'pu69wh78_101117_114745.3.cbin'};
notes_abh_lo = [36 83 36];

% --- jbh
replfiles_jbh_hi = {'pu69wh78_101117_212709.124.cbin', ...
    'pu69wh78_101117_222646.203.cbin', 'pu69wh78_101117_223206.290.cbin'};
notes_jbh_hi = [50 27 64];

replfiles_jbh_lo = {'pu69wh78_101117_104946.142.cbin', ...
    'pu69wh78_101117_114457.188.cbin', 'pu69wh78_101117_114745.3.cbin'};
notes_jbh_lo = [27 51 71];


% ReplacementFiles = {...
%     '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin', ...
%     '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin', ...
%     '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin'}; % list of files - extract replacement syls (provide in order, of sylstoremove)
% NotesToTake = [16 32 47]; % array, corresponding to Replacement files



% ========================== 1) ab(h) high; jb(h) low
ReplacementFiles = [replfiles_abh_hi replfiles_jbh_lo];
NotesToTake = [notes_abh_hi notes_jbh_lo];
suffix = 'aHI_jLO'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)



% ========================== 1) ab(h) high; jb(h) high
ReplacementFiles = [fliplr(replfiles_abh_hi) replfiles_jbh_hi];
NotesToTake = [fliplr(notes_abh_hi) notes_jbh_hi];
suffix = 'aHI_jHI'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)


% ========================== 1) ab(h) lo; jb(h) high
ReplacementFiles = [replfiles_abh_lo fliplr(replfiles_jbh_hi)];
NotesToTake = [notes_abh_lo fliplr(notes_jbh_hi)];
suffix = 'aLO_jHI'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)



% ========================== 1) ab(h) lo; jb(h) lo
ReplacementFiles = [fliplr(replfiles_abh_lo) fliplr(replfiles_jbh_lo)];
NotesToTake = [fliplr(notes_abh_lo) fliplr(notes_jbh_lo)];
suffix = 'aLO_jLO'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)










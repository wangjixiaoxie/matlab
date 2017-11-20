function DATSTRUCT = lt_neural_BatchSmthEXTRACT(basedir, subdirs, Batchfiles, chanstoplot, motifstoplot, ...
    sylstoalign, pretime, posttime, plotRaw, Conditions)
%% lt 11/14/17 - extracts and saves smoothed firing rate for given motifs and channels, for directories and batches

% % ----- different datasets, each with different batch/dir combination
% basedir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111217_Morning_DirUndir/';
% subdirs = {'UNDIR', 'DIR'}; % these must be dirs within basedir. they will be used as fieldnames in extracted structure
% Batchfiles = {'batchall', 'batchall'}; % must be one for each subdir;
% 
% % ------ params for all subdirs
% chanstoplot = [9 14 21]; % chip channels. leave empty to get all
% motifstoplot = {'abh', 'jbh', 'abhh', 'jbhh'}; % cell arary of motifs, strings
% sylstoalign = [3 3 4 4]; % which syl in that motif? (one for each motif
% pretime = 0.1; % sec, from onset
% posttime = 0.15; % sec, from offset
% plotRaw = 0; % plot each extracted trial (raw, smoothed, and spec)


%% 
if ~exist('Conditions', 'var') % this is for naming fields in struct
    Conditions = subdirs;
end

%% =============== run
cd(basedir);

for i=1:length(subdirs)
    
    songdir = [basedir subdirs{i}];
    % where the songs and batch are located
    batchf = Batchfiles{i}; % contains .rhd names
    cond = Conditions{i};
    
    DATSTRUCT.(cond) = lt_neural_BatchSmth(songdir, batchf, chanstoplot, motifstoplot, ...
    sylstoalign, pretime, posttime, plotRaw);
end

%% ============== save

tstamp = lt_get_timestamp(0);
fname = [basedir 'DATSTRUCT_BatchSm_' tstamp '.mat'];
save(fname, 'DATSTRUCT');
daveclust=load('summary_dict.mat');
daveclust.matlab_labels=squeeze(daveclust.matlab_labels);

%% plot all song spectrograms - plot matlab labels and daveclust labels
% Plot cbin same way we do for evsonganaly

% ==== get list of song files daveclust used
songnos=unique(daveclust.syllable_song_number);
songnos=songnos+1;

% === initiate figure
count=1;
SubplotsPerFig=3;
subplotrows=3;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];


% === for each song, load, plot spectrogram, and overlay labels
cc=0;
for i=1:length(songnos);
    ind=songnos(i);
    songname_wav=daveclust.song_files(i,:);
    songname_cbin=[songname_wav(1:end-4) '.cbin'];
    songname_notmat=[songname_cbin '.not.mat'];
    
    % ==== PLOT
    [dat, Fs, DOFILT, ext]=ReadDataFile(songname_cbin,'0');
    [sm,sp,t,f]=SmoothData(dat,Fs,DOFILT,'‘hanningfirff’'); % ends up doing buttor, with filtfilt. (because there is extra quote marks in filter type)
    sp=sp(10:140,:);
    f=f(10:140);
    
    mntmp = min(min(sp(sp>0)));
    sp(sp==0) = mntmp;

    sptemp=log(sp);
    sptemp = sptemp - min(min(sptemp));

    % Plot
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(songname_cbin);
    imagesc(t, f, sptemp);
%     hx=get(gca, 'XLabel'); % to moxe xlabel down, not working yet
%     set(hx, 'Position', get(hx, 'Position')-[100 6000 100])
    axis([t(1) t(end) f(1) f(end)]);
    

    % ==== LABELS (notmat)
    load(songname_notmat);
    
    for j=1:length(onsets)
        lt_plot_text(onsets(j)/1000, -900, labels(j),'r')
    end
    cc=cc+length(onsets);
    
    % ==== LABELS (daveclust)
    inds_this_song=daveclust.syllable_song_number==ind-1;
    onsets=daveclust.onsets(inds_this_song);
    labels=daveclust.labels(inds_this_song);
    
    for j=1:length(onsets)
        lt_plot_text(onsets(j)/1000, -1400, int2str(labels(j)),'b')
    end
    
end

%% === plot individual syllables + spectrograms

% === initiate figure
count=1;
SubplotsPerFig=15;
subplotrows=3;
subplotcols=5;
fignums_alreadyused=[];
hfigs=[];

% params
numsyls=50;
Fs=32000;
DOFILT=1;

for i=1:numsyls
    dat=sylables_raw{i};
    
    [sm,sp,t,f]=SmoothData(dat,Fs,DOFILT,'‘hanningfirff’'); % ends up doing buttor, with filtfilt. (because there is extra quote marks in filter type)
    sp=sp(10:140,:);
    f=f(10:140);
    
    mntmp = min(min(sp(sp>0)));
    sp(sp==0) = mntmp;

    sptemp=log(sp);
    sptemp = sptemp - min(min(sptemp));

    % Plot
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(songname_cbin);
    imagesc(t, f, sptemp);
%     hx=get(gca, 'XLabel'); % to moxe xlabel down, not working yet
%     set(hx, 'Position', get(hx, 'Position')-[100 6000 100])
    axis([t(1) t(end) f(1) f(end)]);
end
   

%% === plot individual syllables + psds

% === initiate figure
count=1;
SubplotsPerFig=15;n
subplotrows=3;
subplotcols=5;
fignums_alreadyused=[];
hfigs=[];

% params
numsyls=50;
Fs=32000;
DOFILT=1;

for i=1:numsyls
    dat=segedpsds(i,:);

    % Plot
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    title(songname_cbin);

    plot(1:length(dat), dat, '-')
    
end
   


%% figuring out best amount to pad beginning and end of syllable




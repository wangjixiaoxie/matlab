%% functions to analyze/extract smoothed firing (as opposed to sorted units)
%% lt 11/10/17 -
clear all; close all;
DATSTRUCT = struct;

% ============== 1) UNDIR
songdir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/MorningDIRtest/UNDIR';
% where the songs and batch are located
batchf = 'batchall'; % contains .rhd names
chanstoplot = [9 14 21]; % chip channels. leave empty to get all
motifstoplot = {'abh', 'jbh', 'abhh', 'jbhh'}; % cell arary of motifs, strings
sylstoalign = [3 3 4 4]; % which syl in that motif? (one for each motif
pretime = 0.1; % sec, from onset
posttime = 0.15; % sec, from offset
plotRaw = 0; % plot each extracted trial (raw, smoothed, and spec)

DATSTRUCT.UNDIR = lt_neural_BatchSmth(songdir, batchf, chanstoplot, motifstoplot, ...
    sylstoalign, pretime, posttime, plotRaw);


% ================= DIR
songdir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/MorningDIRtest/DIR';
% where the songs and batch are located
batchf = 'batchall'; % contains .rhd names
chanstoplot = [9 14 21]; % chip channels. leave empty to get all
motifstoplot = {'abh', 'jbh', 'abhh', 'jbhh'}; % cell arary of motifs, strings
sylstoalign = [3 3 4 4]; % which syl in that motif? (one for each motif
pretime = 0.1; % sec, from onset
posttime = 0.15; % sec, from offset
plotRaw = 0; % plot each extracted trial (raw, smoothed, and spec)

DATSTRUCT.DIR = lt_neural_BatchSmth(songdir, batchf, chanstoplot, motifstoplot, ...
    sylstoalign, pretime, posttime, plotRaw);



%% ================= plot outputs

figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% ==== three subplots for each motif, dir, undir, and overlay dir and undir
nummotifs = length(DATSTRUCT.UNDIR.motifnum);
numchannels = length(chanstoplot);
for i=1:nummotifs
    motifname = DATSTRUCT.UNDIR.motifnum(i).motifname;
    
    %     = ~cellfun('isempty', DATSTRUCT_DIR.motifnum(i).DatAll)
    
    for cc=chanstoplot
        hsplots = [];
        % ----- UNDIR
        Ymean_UNDIR = [];
        Ysem_UNDIR = [];
        dirfield = 'UNDIR';
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['[' dirfield ']' motifname '-ch' num2str(cc)]);
        hsplots = [hsplots hsplot];
        datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{cc};
        t = DATSTRUCT.(dirfield).motifnum(i).t;
        % -- plot raw
        plot(t, datmat, 'Color', [0.7 0.7 0.7]);
        % -- overlay median, with 75th tiles
        datmedian = median(datmat,1);
        datCI = prctile(datmat, [75 25]);
        datSEM = lt_sem(datmat);
        shadedErrorBar(t, datmedian, datSEM, {'Color','r'}, 1);
        line([pretime pretime],ylim);
        lt_plot_zeroline;
        
        Ymean_UNDIR = datmedian;
        Ysem_UNDIR = datSEM;
        
        % ---- DIR
        Ymean_DIR = [];
        Ysem_DIR = [];
        dirfield = 'DIR';
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[' dirfield ']' motifname '-ch' num2str(cc)]);
        datmat = DATSTRUCT.(dirfield).motifnum(i).DatAll{cc};
        t = DATSTRUCT.(dirfield).motifnum(i).t;
        % -- plot raw
        plot(t, datmat, 'Color', [0.7 0.7 0.7]);
        % -- overlay median, with 75th tiles
        datmedian = median(datmat,1);
        datCI = prctile(datmat, [75 25]);
        datSEM = lt_sem(datmat);
        shadedErrorBar(t, datmedian, datSEM, {'Color','r'}, 1);
        line([pretime pretime],ylim);
        lt_plot_zeroline;
        
        Ymean_DIR = datmedian;
        Ysem_DIR = datSEM;
        
        % ----- COMBINED
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[DIR(r) AND UNDIR(k)]' motifname '-ch' num2str(cc)]);
        
        shadedErrorBar(t, Ymean_UNDIR, Ysem_UNDIR, {'Color','k'}, 1);
        shadedErrorBar(t, Ymean_DIR, Ysem_DIR, {'Color','r'}, 1);
        line([pretime pretime],ylim);
        lt_plot_zeroline;
        
        % ---
        linkaxes(hsplots, 'xy');
        
    end
end


%% 1) plot extracted motif for batch of songs that are labeled
close all;
batchf = 'batchall';
chanstoplot = [9 14 21];
motiftoplot = 'jbh';
syltoalign = 3; % which syl in that motif?
pretime = 0.1; % sec, from onset
posttime = 0.2; % sec, from offset


% --- extract data for each batch file
fid = fopen(batchf);
fname = fgetl(fid);


figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

fs = 30000; % code makes sure is really 30k. is ok.

% --------- gaussian for smothing
windowsize=0.02; % from -2sd to +2sd
sigma=(windowsize/4)*fs; %
numsamps=4*sigma; % (get 2 std on each side)
alpha= numsamps/(2*sigma); % N/2sigma
gaussFilter = gausswin(numsamps, alpha);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.


DatAll = cell(1,32); % one for each chan, each cell

% ---- pre/post time in samps
pretime = round(pretime*fs);
posttime = round(posttime*fs);

while ischar(fname)
    
    % ----------- 1) load labels
    try
        notdat = load([fname '.not.mat']);
    catch err
        disp('continueing')
        fname = fgetl(fid);
        continue
    end
    
    % --- find all instances of the motif
    motifpos = strfind(notdat.labels, motiftoplot);
    if isempty(motifpos)
        disp('continuing')
        fname = fgetl(fid);
        continue
    end
    
    % ------- 2) load neural data
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fname);
    assert(fs == frequency_parameters.amplifier_sample_rate, 'not 30k');
    disp('loaded');
    for i=1:length(motifpos)
        pos = motifpos(i)+syltoalign-1;
        
        % ---- extract neural data
        onset = notdat.onsets(pos);
        onset = round(fs*onset/1000); % convert to samps
        
        onsetextract = onset -pretime; % pretime and posttime
        offsetextract = onset + posttime;
        
        %         % -- convert to samples
        %         onsetextract = round(onsetextract*fs);
        %         offsetextract = round(offsetextract*fs);
        %
        % ------- extract neural data
        for cc = 1:length(chanstoplot)
            disp(cc);
            chan = chanstoplot(cc);
            dat = amplifier_data([amplifier_channels.chip_channel] == chan, :);
            
            % --- segment
            dat = dat(onsetextract:offsetextract);
            
            % --- filter, rectify, and smooth
            dat = lt_neural_filter(dat, frequency_parameters);
            datfilt = dat;
            
            dat = abs(dat);
            dat = conv(dat, gaussFilter);
            length(dat), length(datfilt)
            % -- clip off edges
            dat = dat(numsamps/2:end-numsamps/2);
            
            if (1)
                % -- plot raw filtered dat
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([fname '-pos' num2str(pos) '-ch' num2str(chan)]);
                t = [1:length(datfilt)]./fs;
                plot(t, datfilt);
                t = [1:length(dat)]./fs;
                plot(t, dat,'-k');
            end
            % =================== add to output
            DatAll{chan} = [DatAll{chan}; dat];
        end
        
        
        if (1)
            % debugging, make sure songs are aligned
            % --- extract dat
            songdat = board_adc_data(1,:);
            datextract = songdat(onsetextract:offsetextract);
            
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([fname '-pos' num2str(pos)]);
            lt_plot_spectrogram(datextract, fs, 0, 0);
            line([pretime pretime], ylim)
            
        end
        
        
        
    end
    fname = fgetl(fid);
    disp(fname);
end


% ============================ FOR EACH CHAN, PLOT NEURAL
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for cc =1:length(chanstoplot)
    chan = chanstoplot(cc);
    
    % --- extract dat
    datmat = DatAll{chan};
    t = [1:size(datmat,2)]./fs;
    
    % --- plot all trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    plot(t, datmat, 'Color', [0.7 0.7 0.7]);
    
    % --- plot mean,sem
    datmean = mean(datmat,1);
    datsem = lt_sem(datmat);
    shadedErrorBar(t, datmean, datsem, {'Color','k'},1);
    
    line([pretime/fs pretime/fs], ylim);
    
    title([motiftoplot '-ch' num2str(chan)]);
    xlim([0 (pretime+posttime)/fs]);
    
end

function DATSTRUCT = lt_neural_BatchSmth(songdir, batchf, chanstoplot, Motifstoplot, ...
    Sylstoalign, pretime, posttime, plotRaw, plotSpecifics)

clear Sylstoalign; % don't need this, using paranthese to define tokens now.

%% lt 11/11/17

% given batch, motif, and desired channels, pulls out aligned,
% smoothed/rectified neural for all trials

%% params
% songdir = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/MorningDIRtest/DIR';
% % where the songs and batch are located
% batchf = 'batchall'; % contains .rhd names
% chanstoplot = [9 14 21]; % chip channels. leave empty to get all
% Motifstoplot = {'abh', 'jbh'}; % cell arary of motifs, strings
% Sylstoalign = [3 3]; % which syl in that motif? (one for each motif
% pretime = 0.1; % sec, from onset
% posttime = 0.2; % sec, from offset
% plotRaw = 0; % plot each extracted trial (raw, smoothed, and spec)


%% initiate

fs = 30000; % code will makes sure is really 30k. is ok.

% --------- gaussian for smothing
windowsize=0.01; % from -2sd to +2sd
sigma=(windowsize/4)*fs; %
numsamps=4*sigma; % (get 2 std on each side)
alpha= numsamps/(2*sigma); % N/2sigma
gaussFilter = gausswin(numsamps, alpha);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.


% ---- pre/post time in samps
pretime = round(pretime*fs);
posttime = round(posttime*fs);


%% RUN

DATSTRUCT = struct;
DATSTRUCT.params.batchf = batchf;
DATSTRUCT.params.songdir = songdir;
DATSTRUCT.params.Motifstoplot = Motifstoplot;
DATSTRUCT.params.pretime = pretime;
DATSTRUCT.params.posttime = posttime;
% DATSTRUCT.params.Sylstoalign = Sylstoalign;
DATSTRUCT.params.chanstoplot = chanstoplot;

%% === ITERATE OVER MOTIFS
for mm = 1:length(Motifstoplot)
    
    % === params
    motiftoplot = Motifstoplot{mm};
    %     syltoalign = Sylstoalign(mm);
    
    % =================== things to collect
    DatAll = cell(1,32); % one for each chan, each cell, OOUTPUT
    DatAllFilt = cell(1,32);
    FFall = [];
    FnamesAll = {};
    posTokenWithinSong = [];
    NoiseTrials = cell(1,32);
    
    
    % === run
    figcount=1;
    subplotrows=5;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    cd(songdir)
    
    % --- extract data for batch file
    fid = fopen(batchf);
    fname = fgetl(fid);
    
    
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
        [tokenExtents, startinds, endinds, matchlabs] = lt_batchsong_regexp(notdat.labels, motiftoplot);           %%
        if isempty(tokenExtents)
            fname = fgetl(fid);
            continue
        end
        assert(size(tokenExtents,1)==1, 'needs to be horizontal');
        
        if isempty(tokenExtents)
            disp('continuing')
            fname = fgetl(fid);
            continue
        end
        
        %             motifpos = strfind(notdat.labels, motiftoplot);
        %         if isempty(motifpos)
        %             disp('continuing')
        %             fname = fgetl(fid);
        %             continue
        %         end
        %
        % ------- 2) load neural data
        [amplifier_data,board_dig_in_data,frequency_parameters, ...
            board_adc_data, board_adc_channels, amplifier_channels, ...
            board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fname);
        assert(fs == frequency_parameters.amplifier_sample_rate, 'not 30k');
        disp('loaded');
        
        % ---------- 3) try to load file indicating what syl rends are
        % noisy
        if exist([fname '.noise.mat'], 'file')
            trialsToRemove = load([fname '.noise.mat']);
            % --- confirm that this matches current note dat
            if ~all(trialsToRemove.trialsToRemove.notdat.labels ==  notdat.labels)
                disp('PROBLEM - trialsToRemove does not match note dat..., not using ...');
                trialsToRemove = [];
            end
        else
        trialsToRemove = [];    
        end
        
        
        for i=1:length(tokenExtents)
            %             pos = motifpos(i)+syltoalign-1;
            pos = tokenExtents(i);
            % ======================= extract FF
            if exist([fname '.calcff.mat'], 'file')
                load([fname '.calcff.mat']);
                ff = FFstruct.FFall(pos);
            else
                ff = nan;
            end
            FFall = [FFall; ff];
            
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
            t= [];
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
                % -- clip off edges
                dat = dat(numsamps/2:end-numsamps/2);
                t = [1:length(dat)]./fs;
                
                if plotRaw==1
                    
                    % -- only plot specifics
                    if exist('plotSpecifics', 'var')
                        if strcmp(plotSpecifics{1}, motiftoplot) ...
                                & plotSpecifics{2}==chan
                            
                            % -- plot raw filtered dat
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            title([fname '-pos' num2str(pos) '-ch' num2str(chan) '-' motiftoplot]);
                            t = [1:length(datfilt)]./fs;
                            plot(t, datfilt);
                            t = [1:length(dat)]./fs;
                            plot(t, dat,'-k', 'LineWidth', 2);
                            axis tight
                            line([pretime pretime]./fs, ylim);
                        end
                    end
                    
                end
                
                % =================== add to output
                DatAll{chan} = [DatAll{chan}; single(dat)];
                DatAllFilt{chan} = [DatAllFilt{chan}; single(datfilt)];
                
                % ================== note down whether this is a noise file
                if any(ismember(trialsToRemove.trialsToRemove.chan{chan}, pos))
                    % then this is noise trial
                    NoiseTrials{chan} = [NoiseTrials{chan}; 1];
                else
                    % then is not noise
                    NoiseTrials{chan} = [NoiseTrials{chan}; 0];
                end
            end
            
            
            % -- store filename, position in file
            FnamesAll = [FnamesAll fname];
            posTokenWithinSong = [posTokenWithinSong pos];
            
            if plotRaw==1
                % debugging, make sure songs are aligned
                % --- extract dat
                songdat = board_adc_data(1,:);
                datextract = songdat(onsetextract:offsetextract);
                
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([fname '-pos' num2str(pos) '-' motiftoplot]);
                lt_plot_spectrogram(datextract, fs, 0, 0);
                line([pretime pretime]./fs, ylim)
                
            end
            
            
            
        end
        fname = fgetl(fid);
        disp(fname);
    end
    
    % ================ SAVE TO OUTPUT STRUCT
    
    DATSTRUCT.motifnum(mm).motifname = motiftoplot;
    %     DATSTRUCT.motifnum(mm).syltoalign = syltoalign;
    DATSTRUCT.motifnum(mm).DatAll = DatAll;
    DATSTRUCT.motifnum(mm).DatAllRaw = DatAllFilt;
    DATSTRUCT.motifnum(mm).t = single(t);
    DATSTRUCT.motifnum(mm).FF = FFall;
    DATSTRUCT.motifnum(mm).fs = fs;
    DATSTRUCT.motifnum(mm).pretime = pretime;
    
    DATSTRUCT.motifnum(mm).SongnamesAll = FnamesAll;
    DATSTRUCT.motifnum(mm).posTokenWithinSong = posTokenWithinSong;
    DATSTRUCT.motifnum(mm).NoiseTrials = NoiseTrials;
    
end



%% ============================ FOR EACH CHAN, PLOT NEURAL
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
    if size(datmat,1)>1
        shadedErrorBar(t, datmean, datsem, {'Color','k'},1);
    end
    
    line([pretime/fs pretime/fs], ylim);
    
    title([motiftoplot '-ch' num2str(chan)]);
    xlim([0 (pretime+posttime)/fs]);
    
end

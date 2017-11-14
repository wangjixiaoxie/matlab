function lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)

%% defualts
rolltime = 0.002; % s
sidepad = 0.015; % durign splicing
normAmpl = 1;


%% lt 11/10/17 - given a song, splices in syllables at specified locations
% uses evsonganalys segmentation. matinains overall tempo of original file.
% aligns all splices by onsets. if segment is longer then overlaps at ofset. if
% shorter then pads segment with 0.

% saves wav file, with L channel as song, and R channel as digital pulses (of the new song).

% ========= 1) load original song and find locations for splicing out
[dat, fs] = ReadCbinFile(songfname);
datOriginal = dat;


notdat = load([songfname '.not.mat']);

% ================== iterate splice thru each location
for i=1:length(sylstoremove)
    
    % =========== 1) EXTRACT SONG
    sylnum = sylstoremove(i);
    sampOn = round(fs*(notdat.onsets(sylnum)/1000 - sidepad));
    sampOff = round(fs*(notdat.offsets(sylnum)./1000 + sidepad));
    
    
    % =========== 2) extract replacement syl, and do roll on/off
    [datReplace, fs2] = ReadCbinFile(ReplacementFiles{i});
    
%     % --- bandpass filter
% filter_type='hanningfirff';
% F_low  = 500;
% F_high = 8000;
% datReplace=bandpass(datReplace,fs,F_low,F_high,filter_type);

    
    %     datReplace = datReplace/max(datReplace);
    assert(fs2==fs, 'asdfas');
    notdatReplace = load([ReplacementFiles{i} '.not.mat']);
    
    % -- ons and off
    sylnumRepl = NotesToTake(i);
    onRepl = round(fs*(notdatReplace.onsets(sylnumRepl)/1000 - sidepad));
    offRepl = round(fs*(notdatReplace.offsets(sylnumRepl)/1000 + sidepad));
    
    segRepl = datReplace(onRepl:offRepl);
    
    % ============ normalize amplitude if desired
    if normAmpl ==1
    segampl = prctile(log(segRepl.^2  + 0.0001), 75);
    datampl = prctile(log(dat(sampOn:sampOff).^2 + 0.0001), 75);
    
    segRepl = segRepl*(datampl/segampl);
    end
    
    % ============= 3) make size of segment match size of thing taking out.
    % --- if size not the same, align by onset. if segment is shorter, add
    % zeros.
    tmp = (sampOff - sampOn+1) - length(segRepl);
    
    % if this is large, or larger than sidepad, then error
    if abs(tmp/min([(sampOff - sampOn+1), length(segRepl)])) > 0.1
        disp('PROBLEM, >10% duration diff')
    end
    if (tmp/fs)>sidepad
        disp('RPBLEM, dur diff graeter than side pad')
    end
    
    if tmp>0
        % then seg too short, pad with zeros
        segRepl = [segRepl; zeros(tmp,1)];
    elseif tmp<0
        % then too long, make sampOff larger
        sampOff = sampOff-tmp;
    end
    
    
    
    % =========== 4)  do all roll offs/ons
    expsize = round(rolltime*fs);
    
    % -------- 1) DATA
    % --- roll off (left of edge)
    smth = exp(-((1:expsize)-1)/(expsize/5));
    dat(sampOn-expsize:sampOn-1) = dat(sampOn-expsize:sampOn-1).*smth';
    
    % -- roll on (right of edge)
    smth = exp(-((expsize:-1:1)-1)/(expsize/5));
    dat(sampOff+1:sampOff+expsize) = dat(sampOff+1:sampOff+expsize).*smth';
    
    % ---------------- 2) SEGMMENT TO INPUT
    % --- roll on
    smth = exp(-((expsize:-1:1)-1)/(expsize/5));
    segRepl(1:expsize) = segRepl(1:expsize).*smth';
    
    % --- roll off
    smth = exp(-((1:expsize)-1)/(expsize/5));
    segRepl(end-expsize+1:end) = segRepl(end-expsize+1:end).*smth';
    
    if (0)
        figure; hold on;
        subplot(211);
        plot(segRepl);
        
        subplot(212); hold on;
        lt_plot_spectrogram(segRepl, fs, 1, 0);
    end
    
    
    
    % ============= insert segment into song
    dat(sampOn:sampOff) = segRepl;
    
    if (1)
        lt_figure; hold on;
        subplot(411);
        plot(datOriginal, 'k');
        line([sampOn sampOn], ylim);
        
        line([sampOff sampOff], ylim);
        
        subplot(412);
        lt_plot_spectrogram(datOriginal, fs, 1, 0);
        set(gca, 'Ydir', 'reverse')
        
        subplot(413); hold on;
        plot(dat, 'c');
        plot(sampOn:sampOff, segRepl, 'r');
        
        subplot(414); hold on;
        
        lt_plot_spectrogram(dat, fs, 1, 0);
    end
end


% ============= extract onsets and offsets to save as well
% for now, use the onsets and offsets for the original file. the onsets
% should not change, and the offsets should change only slightly
datdig = zeros(length(dat),1);

if UseOriginalDigPulse==1
    onsets = notdat.onsets;
    offsets = notdat.offsets;
    
else
    % then makes the new digital pulses for the new song, using same
    % segmenetation algorithm
    sm=SmoothData(dat,fs,1,'hanningfirff');
    sm(1)=0.0;sm(end)=0.0;
    [ons,offs]=SegmentNotesJC(sm,fs,notdat.min_int,notdat.min_dur,notdat.threshold);
    onsets=ons*1e3;offsets=offs*1e3;
end

for i=1:length(onsets)
    ons = onsets(i);
    off = offsets(i);
    
    % convert from ms to samps
    ons = round(fs*(ons/1000));
    off = round(fs*(off/1000));
    
    datdig(ons:off) = 1;
end


% ================ filter and scale by max
% --- bandpass filter
filter_type='hanningfirff';
F_low  = 500;
F_high = 8000;
dat=bandpass(dat,fs,F_low,F_high,filter_type);

% -- scale
dat = dat./max(dat);

%% save
% ############# 1) forward - combine song and dig signal
tmp = strfind(songfname, '/');
fnameout = songfname(tmp(end)+1:end);
fnameout = [fnameout '_' suffix '.wav'];
fnameout = [OutDir '/' fnameout];

wavwrite([dat datdig], fs, 16, fnameout);

% ############# MAKE REVERSE VERSION
dat_rev = flipud(dat);
datdig_rev = flipud(datdig);

tmp = strfind(songfname, '/');
fnameout = songfname(tmp(end)+1:end);
fnameout = [fnameout '_' suffix '_REV.wav'];
fnameout = [OutDir '/' fnameout];
wavwrite([dat_rev datdig_rev], fs, 16, fnameout);



% ########### plot
if (1)
    figure; hold on;
    subplot(211); hold on;
    title('forward');
    plot(datOriginal./max(datOriginal), 'k');
    plot(dat, 'c');
    plot(datdig, 'r');
    
    subplot(212); hold on;
    title('rev');
    plot(dat_rev, 'c');
    plot(datdig_rev, 'r');
end

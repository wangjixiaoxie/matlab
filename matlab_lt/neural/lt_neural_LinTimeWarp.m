%% - LT, linearly time warp motifs, project neural to that
% assumes that onset of syl 1 to offset of syl 2 should be aligned.

function [SegmentsExtract, Params] = lt_neural_LinTimeWarp(SegmentsExtract, Params, TargetDur)

% SongDat
% SegmentsExtract
% Params
% all those params as outputs from lt_neural_RegExp
% TargetDur = []; then uses median dur across all trials;
% TargetDur = 0.5s; then forces 0.5s (onset(syl1) to offset(syllast) (i.e.
% ignoreing pre and post dur);


numtrials=length(SegmentsExtract);
predur=Params.REGEXP.predur;
postdur=Params.REGEXP.postdur;

%% ---- SIMPLE - linearly stretch entire motif (IGNORES predur and postdur)

% === 1) find trial of median length
motifontimes = predur + ...
    [SegmentsExtract.global_ontime_motifInclFlank];

motifofftimes= [SegmentsExtract.global_offtime_motifInclFlank] - ...
    postdur;

allmotiftimes=motifofftimes-motifontimes;

if isempty(TargetDur)
    % then find median duration and use that
    motiftime_med=median(allmotiftimes);
else
    % use
    motiftime_med= TargetDur;
end

% - plot distrib
lt_figure; hold on; lt_plot_histogram(allmotiftimes)
line([motiftime_med motiftime_med], ylim);
title('distribution of motif durs (exclude flank)');
xlabel('sec');

% =========== 2) for each trial, stretch or compress spike times so that
% aligned with median.
% note: all trials will still be aligned at time=predur.
lt_figure; hold on;
title('bk-old; rd - scaled');
ylabel('trial'); xlabel('sec');

for i=1:numtrials
    disp(['rescaling spikes trial ... ' num2str(i)]);
    % - get motif duration
    motifon = predur + SegmentsExtract(i).global_ontime_motifInclFlank;
    motifoff = SegmentsExtract(i).global_offtime_motifInclFlank - postdur;
    
    motifdur=motifoff-motifon;
    %     plot(motifdur, 1, 'o');
    
    scaleFactor=motifdur/motiftime_med; % >1 means need to compress spiketimes
    
    % --- scale spike times.
    spktimes=SegmentsExtract(i).spk_Times;
    spktimes_rezero=spktimes - predur; % subtract constant to make the aligned time (predur) equal to 0
    spktimes_rezero_scaled = spktimes_rezero./scaleFactor;
    spktimes_scaled=spktimes_rezero_scaled + predur; % add on the predur time, so that they are now all aligned at predur again.
    
    % save
    SegmentsExtract(i).spk_Times_scaled=spktimes_scaled;
    SegmentsExtract(i).scaleFactor=scaleFactor;
    
    % - plot spike times
    numspks=length(SegmentsExtract(i).spk_Times);
    
    for j=1:numspks
        spktime=SegmentsExtract(i).spk_Times(j);
        line([spktime spktime], [i-0.3 i+0.3], 'Color', 'k');
        
        spktime2=SegmentsExtract(i).spk_Times_scaled(j);
        line([spktime2 spktime2], [i-0.3 i+0.3], 'Color', 'r');
    end
    
    lt_plot_text(0, i+0.2, ['scaleFactor=' num2str(scaleFactor)], 'g');
    
    % ===== SCALE SONG DATA (interpolate)
    % IN PROGRESS - would need to determine new inds, then interpolate
    % those inds (with aligned point at 0)
%     songdat=SegmentsExtract.songdat;
%     songdat_rezero
    
    
end
line([predur predur], ylim, 'Color', 'g');





%% 
Params.LINTIMEWARP.motiftime_median=motiftime_med;



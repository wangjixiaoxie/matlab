function SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum, ...
    kernelSD, binsize_spks, extractSingleTrials, segextract_for_trialdur)

%%
% in case the seg extract for dat doesnt have timing variables, enter
% another one that does ...
if ~exist('segextract_for_trialdur', 'var')
    segextract_for_trialdur = [];
else
    assert(length(segextract_for_trialdur) == length(SegmentsExtract), 'sdaf');
end


%%

TryDiffMethods = 0; % if 1, then overlays result sof diff FR estimation methods. if 0, then uses gaussian kernel smoothing
% note: will automaticalyl plot figure as well.

% if want all clusters, then set clustnum = [];

if ~exist('extractSingleTrials', 'var')
    extractSingleTrials=0; % if 0, then will only extract one large FR mat with common duration. If 1, then cell, with each potentiayl different length
elseif isempty(extractSingleTrials)
    extractSingleTrials = 0;
end

%% input SegmentsExtract. output SegmentsExtract, but with smoothed FR as a field
% uses kde, gaussian.

if ~exist('kernelSD', 'var')
    kernelSD = 0.005; % ms gaussian
elseif isempty(kernelSD)
    kernelSD = 0.005; % ms gaussian
end

if ~exist('binsize_spks', 'var')
    binsize_spks = 0.001; % for binning spikes
elseif isempty(binsize_spks)
    binsize_spks = 0.001;
end

%%
% clustnum = desired cluster
if isfield(SegmentsExtract, 'spk_Times') % only try to smooth if there is any data
    
    NumTrials = length(SegmentsExtract);
    
    if ~isempty(segextract_for_trialdur)
        commonTrialDur = min([segextract_for_trialdur.global_offtime_motifInclFlank] - ...
            [segextract_for_trialdur.global_ontime_motifInclFlank]); % shortest trial
    else
        commonTrialDur = min([SegmentsExtract.global_offtime_motifInclFlank] - ...
            [SegmentsExtract.global_ontime_motifInclFlank]); % shortest trial
    end
    maxInd = floor(commonTrialDur/binsize_spks);
    
    for i=1:NumTrials
        if isfield(SegmentsExtract, 'spk_Clust')
            % then is my data -, match to clust
            if isempty(clustnum)
                inds = SegmentsExtract(i).spk_Clust>0;
            else
                inds = SegmentsExtract(i).spk_Clust == clustnum;
            end
            spktimes = SegmentsExtract(i).spk_Times(inds);
        else
            % is Sober/Mel data
            spktimes = SegmentsExtract(i).spk_Times;
        end
        
        
        if (0)
            disp('NOTE: KDE NOT YET WRITTEN, USING BOXCAR SMOOTHING');
            
            [xbin, ~, ~, ymeanhz] = lt_neural_plotRastMean({spktimes}, 0.0075, 0.002, 0, '');
            
            SegmentsExtract(i).FRsmooth_rate = ymeanhz;
            SegmentsExtract(i).FRsmooth_xbin = xbin;
            
        else
            %%
            %             disp('KERNEL SMOOTHING')
            if ~isempty(segextract_for_trialdur)
                trialdur = segextract_for_trialdur(i).global_offtime_motifInclFlank - ...
                    segextract_for_trialdur(i).global_ontime_motifInclFlank;
                
            else
                trialdur = SegmentsExtract(i).global_offtime_motifInclFlank - ...
                    SegmentsExtract(i).global_ontime_motifInclFlank;
            end
            
            binedges = 0:binsize_spks:trialdur;
            [bincounts] = histc(spktimes, binedges); % binned spikes
            
            if isempty(spktimes)
                bincounts = zeros(size(binedges));
            end
            
            x = -3*kernelSD:binsize_spks:3*kernelSD;
            kern = (1/(sqrt(2*pi)*kernelSD)) * exp(-x.^2/(2*kernelSD^2)); % kernel
            
            smoothFR = conv(bincounts(1:end-1), kern, 'same'); % output
            
            % --- visualize and compare TO OTHERS
            if TryDiffMethods==1
                lt_figure; hold on;
                title('bk=kern smth');
                plot(spktimes, 1, 'ok');
                plot(binedges(1:end-1) + binsize_spks/2, smoothFR, '-xk');
                
                
                % -- boxcar
                [xbin, ~, ~, ymeanhz] = lt_neural_plotRastMean({spktimes}, kernelSD*4, 0.002, 0, '');
                plot(xbin, ymeanhz, '-b'); %
                
                % -- optimized kernel bandwidth
                [y, t, optw] = sskernel(spktimes, binedges(1:end-1));
                y = y.*(length(spktimes)/trialdur); % convert from numspks to spkrate (APPROXIMATE!!)
                plot(t, y, 'r-');
                
                lt_plot_text(t(1), y(1), ['opt. bwidth=' num2str(optw)], 'r');
                
                % -- optimized, locally adaptive kernel
                [y, t, optw] = ssvkernel(spktimes, binedges(1:end-1));
                y = y.*(length(spktimes)/trialdur); % convert from numspks to spkrate (APPROXIMATE!!)
                plot(t, y, 'm-');
                lt_plot_text(t(1), y(1)+10, ['local/bwidth; median opt. bwidth=' num2str(median(optw))], 'm');
                
                % --- empirical bayes
                [time, rate, regularity] = BayesRR(spktimes);
                plot(time, rate, 'k--');
                lt_plot_annotation(1, 'bayes (dashed black)', 'k');
                
                % --- bayesian adaptive regression splines
                summary = barsP(bincounts, [0 trialdur], 1);
                plot(binedges(1:end-1), 20*(length(spktimes)/trialdur).*summary.mean(1:end-1), '-g');
                lt_plot_annotation(2,['bars'], 'g');
                
            end
            
            % --- save
            if extractSingleTrials==1
                SegmentsExtract(i).FRsmooth_rate = single(smoothFR');
                SegmentsExtract(i).FRsmooth_xbin = single(binedges(1:end-1)');
            end
            SegmentsExtract(i).FRsmooth_rate_CommonTrialDur = single(smoothFR(1:maxInd)');
            SegmentsExtract(i).FRsmooth_xbin_CommonTrialDur = single(binedges(1:maxInd)');
            
            
        end
    end
end

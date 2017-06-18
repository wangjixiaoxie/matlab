function SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum, kernelSD, binsize_spks)
%%

TryDiffMethods = 0; % if 1, then overlays result sof diff FR estimation methods. if 0, then uses gaussian kernel smoothing
% note: will automaticalyl plot figure as well.


%% input SegmentsExtract. output SegmentsExtract, but with smoothed FR as a field
% uses kde, gaussian.

if ~exist('kernelSD', 'var')
    kernelSD = 0.005; % ms gaussian
end

if ~exist('binsize', 'var')
    binsize_spks = 0.001; % for binning spikes
end

% clustnum = desired cluster
if isfield(SegmentsExtract, 'spk_Clust'); % only try to smooth if there is any data
    
    NumTrials = length(SegmentsExtract);
    
    for i=1:NumTrials
        inds = SegmentsExtract(i).spk_Clust == clustnum;
        spktimes = SegmentsExtract(i).spk_Times(inds);
        
        commonTrialDur = min([SegmentsExtract.global_offtime_motifInclFlank] - ...
                [SegmentsExtract.global_ontime_motifInclFlank]); % shortest trial 
        maxInd = floor(commonTrialDur/binsize_spks);
            
        if (0)
            disp('NOTE: KDE NOT YET WRITTEN, USING BOXCAR SMOOTHING');
            
            [xbin, ~, ~, ymeanhz] = lt_neural_plotRastMean({spktimes}, 0.0075, 0.002, 0, '');
            
            SegmentsExtract(i).FRsmooth_rate = ymeanhz;
            SegmentsExtract(i).FRsmooth_xbin = xbin;
            
        else
            %%
%             disp('KERNEL SMOOTHING')
            
            trialdur = SegmentsExtract(i).global_offtime_motifInclFlank - ...
                SegmentsExtract(i).global_ontime_motifInclFlank;
            
            binedges = 0:binsize_spks:trialdur;
            [bincounts] = histc(spktimes, binedges); % binned spikes
            
            if isempty(spktimes)
                bincounts = zeros(size(binedges));
            end
            
            x = -3*kernelSD:binsize_spks:3*kernelSD;
            kern = (1/(sqrt(2*pi)*kernelSD)) * exp(-x.^2/(2*kernelSD^2)); % kernel
            
            smoothFR = conv(bincounts(1:end-1), kern, 'same'); % output
            
            % --- visualize and compare TO OTHERS
            if TryDiffMethods==1;
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
            SegmentsExtract(i).FRsmooth_rate = smoothFR';
            SegmentsExtract(i).FRsmooth_xbin = binedges(1:end-1)';
            SegmentsExtract(i).FRsmooth_rate_CommonTrialDur = smoothFR(1:maxInd)';
            SegmentsExtract(i).FRsmooth_xbin_CommonTrialDur = binedges(1:maxInd)';
            
            
        end
    end
end

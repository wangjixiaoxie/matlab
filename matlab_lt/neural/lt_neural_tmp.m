function lt_neural_tmp(MOTIFSTATS, NeuronDatabase, motif_regexpr_str, WNchangeDateStrings)

NumNeurons=length(NeuronDatabase.neurons);
NumMotifs=length(motif_regexpr_str);

window=0.02;
windshift=0.004;
motif_predur=0.1;
motif_postdur=0.1;

OnlyPlotNoHit=0;
TrialBinSize=15;

spktimefield='spk_Times';

%%
for i=1:NumNeurons
    
    figcount=1;
    subplotrows=6;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
hsplots1=[];
hsplots2=[];


    for m=1:NumMotifs
        
        
        segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        clustnum=NeuronDatabase.neurons(i).clustnum;
        
        % --- IF ONLY PLOT NO HIT, THEN FIRST PULL OUT THOSE TRIALS
        if OnlyPlotNoHit==1
            indsWithNoHits=[segextract.hit_WN]~=1;
            segextract=segextract(indsWithNoHits);
        end

        
        % ===================== collect spikes for all trials
        numtrials=length(segextract);
        Yspks={};
        for j=1:numtrials
            inds=segextract(j).spk_Clust==clustnum;
            spktimes=segextract(j).(spktimefield)(inds);
            Yspks{j}=spktimes;
        end
        
        % ===================== BASELINE mean and std spiking
        datestring=WNchangeDateStrings{1}; % assumes this is WN on.
        dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
        inds=[segextract.song_datenum] < dnum; %
        [xbase, ~, ~, ymean_hz_base, ~, ~, ystd_hz_base] = ...
            lt_neural_plotRastMean(Yspks(inds), window, windshift, 0, '');

        
        % ========================= GET RUNNING D-PRIME OVER TIME, IN
        % SPECIFIC PREMOTOR WINDOW
        premotor_wind=[-0.08 0.02]; % in sec, relative to onset of token in motif.
        DprimeVals=[];
        AbsDprimeVals=[];
        MedianDatenumVals = [];
        MedianDayVals=[];
        FFVals=[];
        AreAllTrialsBaseline = [];
        AreAllTrialsDurWN = [];
        
        
        for j=1:(numtrials-TrialBinSize+1);
            trialbins=j:(j+TrialBinSize-1);
            
            [x, ~, ~, ymean_hz, ~, ~, ystd_hz] =  ...
                lt_neural_plotRastMean(Yspks(trialbins), window, windshift, 0, '');
            
            % -- calc d-prime
            xinds=1:min([length(ymean_hz) length(ymean_hz_base)]);
            numer_tmp=ymean_hz(xinds)-ymean_hz_base(xinds);
            denom_tmp=sqrt(0.5*(ystd_hz_base(xinds).^2 + ystd_hz(xinds).^2));
            dprime=numer_tmp./denom_tmp;
            
            % -- slice dprime within premotor window
            xinds_to_keep = x>=(motif_predur + premotor_wind(1) + window/2) ...
                & x<=(motif_predur + premotor_wind(2) - window/2); % accounting for smoothing window size.
            
            dprime_val = mean(dprime(xinds_to_keep));        
            DprimeVals=[DprimeVals dprime_val];
            
            abs_dprime_val = mean(abs(dprime(xinds_to_keep)));        
          AbsDprimeVals=[AbsDprimeVals abs_dprime_val];
            
            % -- get median time for this bin
            median_datenum=median([segextract(trialbins).song_datenum]);
            firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
            eventtime=datestr(median_datenum, 'ddmmmyyyy-HHMM');
            tmp=lt_convert_EventTimes_to_RelTimes(firstday, {eventtime}); % days from start of expt
            
            MedianDatenumVals = [MedianDatenumVals median_datenum];
            MedianDayVals = [MedianDayVals tmp.FinalValue];
            
            % -- FF
            mean_FF=mean([segextract(trialbins).FF_val]);
            FFVals=[ FFVals mean_FF];
            
            % -- all trials baseline?
            datestring=WNchangeDateStrings{1};
            dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
            if segextract(trialbins(end)).song_datenum < dnum 
                AreAllTrialsBaseline = [AreAllTrialsBaseline 1]; % all trials base
                AreAllTrialsDurWN = [AreAllTrialsDurWN 0];
            elseif segextract(trialbins(1)).song_datenum >= dnum
                AreAllTrialsBaseline = [AreAllTrialsBaseline 0]; 
                AreAllTrialsDurWN = [AreAllTrialsDurWN 1];
            else
                AreAllTrialsBaseline = [AreAllTrialsBaseline 0]; 
                AreAllTrialsDurWN = [AreAllTrialsDurWN 0];
            end                
            
            
        end
        
        
        % ==== PLOT
        if (0) % set as 1 to plot mean dprime, otherwise plots mean abs(dprime)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots1=[hsplots1 hsplot];
title(motif_regexpr_str{m}); ylabel('mean premotor dprime');
plot(MedianDayVals, DprimeVals, 'ok');
plot(MedianDayVals(logical(AreAllTrialsBaseline)), DprimeVals(logical(AreAllTrialsBaseline)), 'ob');
plot(MedianDayVals(logical(AreAllTrialsDurWN)), DprimeVals(logical(AreAllTrialsDurWN)), 'or');
        else
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots1=[hsplots1 hsplot];
title(motif_regexpr_str{m}); ylabel('mean abs(premotor dprime)');
plot(MedianDayVals, AbsDprimeVals, 'ok');
plot(MedianDayVals(logical(AreAllTrialsBaseline)), AbsDprimeVals(logical(AreAllTrialsBaseline)), 'ob');
plot(MedianDayVals(logical(AreAllTrialsDurWN)), AbsDprimeVals(logical(AreAllTrialsDurWN)), 'or');            
        end
        

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots2=[hsplots2 hsplot];
title(motif_regexpr_str{m}); ylabel('ff (hz)');
plot(MedianDayVals, FFVals, 'sk');
plot(MedianDayVals(logical(AreAllTrialsBaseline)), FFVals(logical(AreAllTrialsBaseline)), 'sb');
plot(MedianDayVals(logical(AreAllTrialsDurWN)), FFVals(logical(AreAllTrialsDurWN)), 'sr');

        
    end
    
    linkaxes(hsplots1, 'xy');
    linkaxes(hsplots2, 'xy');
    lt_subtitle(['neuron' num2str(i)]);
end

function lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct, PlotRaw)

%% lt 5/3/17 - effect of bout position on structure and neural (Joanne)

%% params
motif_predur = 0.1;
motif_postdur = 0.02;
PlotSmallTimeWindow = 1;

window = 0.01; % neural;
windshift = 0.002;
window_meanFF = [-0.07 0.02]; % relative onset


%% === for each neuron, extract data for each reg expr string from actual motifs
NumBirds = length(SummaryStruct.birds);
MOTIFSTATS = struct;

for i=1:NumBirds
    NumNeurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:NumNeurons
        
        cd(SummaryStruct.birds(i).neurons(ii).dirname);
        
        %-- load dat
        batchf=SummaryStruct.birds(i).neurons(ii).batchfilename;
        channel_board=SummaryStruct.birds(i).neurons(ii).channel;
        extractSound = 0;
        cd ..
        [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, extractSound);
        
        
        % ---- EXTRACT DAT FOR EACH MOTIF
        motif_regexpr_str = SummaryStruct.birds(i).neurons(ii).POSTINFO.MotifsActual_regexpStr;
        NumSyls = length(motif_regexpr_str);
        
        for j=1:NumSyls
            
            regexpr_str = motif_regexpr_str{j};
            
            alignByOnset=1;
            WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
            % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
            FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
            FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
            % +1 is 1 after token
            FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
            collectWNhit = 0;
            collectWholeBoutPosition=1; % to get position f motif in each bout.
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                regexpr_str, motif_predur, motif_postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
                0,0, collectWNhit, collectWholeBoutPosition);
            
            % ================== OUTPUT
            MOTIFSTATS.birds(i).neurons(ii).motif(j).SegmentsExtract = SegmentsExtract;
            MOTIFSTATS.birds(i).neurons(ii).motif(j).Params = Params;
            MOTIFSTATS.birds(i).neurons(ii).motif(j).motif = regexpr_str;
            
        end
        MOTIFSTATS.birds(i).neurons(ii).clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
    end
end

%% ==== PLOT
spktimefield = 'spk_Times';
% ==== 1) FOR EACH NEURON,
if PlotRaw ==1
    NumBirds = length(SummaryStruct.birds);
    for i=1:NumBirds
        
        
        NumNeurons = length(SummaryStruct.birds(i).neurons);
        rowcount = 1;
        numrows = 4;
        for ii=1:NumNeurons
            
            NumMotifs = length(MOTIFSTATS.birds(i).neurons(ii).motif);
            
            for j=1:NumMotifs;
                
                if rowcount==1;
                    lt_figure; hold on;
                end
                SegmentsExtract = MOTIFSTATS.birds(i).neurons(ii).motif(j).SegmentsExtract;
                if ~isfield(SegmentsExtract, 'FF_val')
                    continue
                end
                FFvals = [SegmentsExtract.FF_val];
                
                if all(isnan(FFvals))
                    continue
                end
                RendInBout = [SegmentsExtract.BOUT_RendInBout];
                
                % FR
                numtrials = length(SegmentsExtract);
                clustnum = MOTIFSTATS.birds(i).neurons(ii).clustnum;
                Yspks={};
                
                for k=1:numtrials
                    inds=SegmentsExtract(k).spk_Clust==clustnum;
                    spktimes=SegmentsExtract(k).(spktimefield)(inds);
                    
                    % look at a small time window
                    windowtmp_forspks = [0 motif_predur+0.1]; % to onset + 0.15;
                    spktimes = spktimes(spktimes>windowtmp_forspks(1) & spktimes<windowtmp_forspks(2));
                    
                    Yspks=[Yspks spktimes];
                end
                % % -- convert to smoothed rate
                % [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
                
                maxrendInbout = max(RendInBout);
                plotcols = lt_make_plot_colors(maxrendInbout, 0, [1 0 0]);
                colnum = 5;
                usedcells = (rowcount-1)*(colnum+2);
                for k=1:maxrendInbout
                    inds = find(RendInBout == k);
                    
                    ffvalstmp = FFvals(inds);
                    yspkstmp = Yspks(inds);
                    
                    % --- get smoothed FR
                    [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(yspkstmp, window, windshift, 0, '');
                    
                    % --- get mean FR in a certain window
                    windowtmp = window_meanFF(2)-window_meanFF(1);
                    windowshifttmp = motif_predur+window_meanFF(1);
                    [xtmp, ~, ~, ytmp, ysemtmp] = lt_neural_plotRastMean(yspkstmp, windowtmp, windowshifttmp, 0, '');
                    FRmean = ytmp(2);
                    FRsem = ysemtmp(2);
                    
                    % ==== PLOT
                    % 1) ffvals
                    lt_subplot(numrows,colnum+2, 1+usedcells); hold on;
                    %     title('FF'); xlabel('pos in bout'); ylabel('mean FF, sem');
                    lt_plot(k, mean(ffvalstmp), {'Errors', lt_sem(ffvalstmp), 'Color', plotcols{k}});
                    
                    % 2) smoothed nerual
                    lt_subplot(numrows, colnum+2, [3:colnum+2]+usedcells); hold on;
                    %     title('neural');
                    X = xbin + (k-1)*windowtmp_forspks(end);
                    shadedErrorBar(X, ymean_hz, ysem_hz, {'Color', plotcols{k}}, 1);
                    xlinepos = motif_predur + (k-1)*windowtmp_forspks(end);
                    line([xlinepos xlinepos], [0 100], 'Color','k');
                    line([xlinepos+window_meanFF(1) xlinepos+window_meanFF(1)], [0 100], 'Color','b', 'LineStyle', '--'); % for mean FF dat
                    line([xlinepos+window_meanFF(2) xlinepos+window_meanFF(2)], [0 100], 'Color','b', 'LineStyle', '--');
                    %     grid on;
                    %     lt_plot_zeroline;
                    %     axis tight
                    
                    % 3) FRmean
                    lt_subplot(numrows, colnum+2, 2+usedcells); hold on;
                    %     title('neural');
                    lt_plot(k, FRmean, {'Errors', FRsem, 'Color', plotcols{k}});
                    %     lt_plot_zeroline;
                end
                
                % minor stuff
                lt_subplot(numrows,colnum+2,1+usedcells); hold on;
                %     title('FF'); xlabel('pos in bout'); ylabel('mean FF, sem');
                axis tight;
                
                lt_subplot(numrows, colnum+2, [3:colnum+2]+usedcells); hold on;
                %     title('neural');
                grid on;
                axis tight
                lt_plot_zeroline;
                lt_plot_annotation(2, ['bird' num2str(i) ', neuron' num2str(ii) ', ' ...
                    MOTIFSTATS.birds(i).neurons(ii).motif(j).motif], 'k');
                
                lt_subplot(numrows, colnum+2, 2+usedcells); hold on;
                %     title('neural');
                axis tight
                
                rowcount = mod(rowcount, numrows)+1;
            end
        end
    end
end

%% ========== COLLECT ALL DATA - PLOT 1) PITCH, 2) FR, 3) CORRELATION
% [MEANS]


BirdnumAll=[];
NeuronNumAll = [];
MotifNumAll = [];
FFmeanAll = [];
FFsemAll = [];
FRmeanAll = [];
FRsemAll = [];
RendInBoutAll = [];
SampSizeAll = [];

NumBirds = length(SummaryStruct.birds);
for i=1:NumBirds
    
    NumNeurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:NumNeurons
        
        NumMotifs = length(MOTIFSTATS.birds(i).neurons(ii).motif);
        
        for j=1:NumMotifs;
            
            SegmentsExtract = MOTIFSTATS.birds(i).neurons(ii).motif(j).SegmentsExtract;
            
            % ==== EXTRACT STUFF
            % -- 1) ff
            if ~isfield(SegmentsExtract, 'FF_val')
                continue
            end
            FFvals = [SegmentsExtract.FF_val];
            
            if all(isnan(FFvals))
                continue
            end
            
            % -- 1b) rend pos
            RendInBout = [SegmentsExtract.BOUT_RendInBout];
            
            % 2) FR
            numtrials = length(SegmentsExtract);
            clustnum = MOTIFSTATS.birds(i).neurons(ii).clustnum;
            Yspks=cell(1, numtrials);
            for k=1:numtrials
                inds=SegmentsExtract(k).spk_Clust==clustnum;
                spktimes=SegmentsExtract(k).(spktimefield)(inds);
                
                % look at a small time window
                windowtmp_forspks = [0 motif_predur+0.1]; % to onset + 0.15;
                spktimes = spktimes(spktimes>windowtmp_forspks(1) & spktimes<windowtmp_forspks(2));
                
                Yspks{k} = spktimes;
                
            end
            
            maxrendInbout = max(RendInBout);
            for k=1:maxrendInbout
                inds = find(RendInBout == k);
                
                if length(inds)<3
                    continue
                end
                
                ffvalstmp = FFvals(inds);
                yspkstmp = Yspks(inds);
                
                % --- get smoothed FR
                if (0)
                    [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(yspkstmp, window, windshift, 0, '');
                end
                
                % --- get mean FR in a certain window
                windowtmp = window_meanFF(2)-window_meanFF(1);
                windowshifttmp = motif_predur+window_meanFF(1);
                [xtmp, ~, ~, ytmp, ysemtmp] = lt_neural_plotRastMean(yspkstmp, windowtmp, windowshifttmp, 0, '');
                
                if length(ytmp)<2
                    continue
                end
                FRmean = ytmp(2);
                FRsem = ysemtmp(2);
                
                % =============== OUTPUT
                BirdnumAll=[BirdnumAll i];
                NeuronNumAll = [NeuronNumAll ii];
                MotifNumAll = [MotifNumAll j];
                FFmeanAll = [FFmeanAll mean(ffvalstmp)];
                FFsemAll = [FFsemAll lt_sem(ffvalstmp)];
                FRmeanAll = [FRmeanAll FRmean];
                FRsemAll = [FRsemAll FRsem];
                assert(all(RendInBout(inds) == k), 'asfasdf');
                RendInBoutAll = [RendInBoutAll k];
                SampSizeAll = [SampSizeAll length(ffvalstmp)];
            end
            
        end
    end
end


%% ====================== PLOT
% ---- 1) for each bird/neuron, overlay all motifs
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
minNumRends = 3; % i.e. at least this many positions in motif.


for i=1:max(BirdnumAll)
    for ii=1:max(NeuronNumAll)
        
        % --- PITCH
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i) ' neuron' num2str(ii)]);
        ylabel('ff (minus base');
        plotcols = lt_make_plot_colors(max(MotifNumAll), 0, 0);
        for j=1:max(MotifNumAll)
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            disp([num2str(RendInBoutAll(inds)) '-' num2str(MotifNumAll(inds))]);
            
            rendsinbout = RendInBoutAll(inds);
            motifnum = MotifNumAll(inds);
            ffmeans = FFmeanAll(inds);
            ffsems = FFsemAll(inds);
            fratemeans = FRmeanAll(inds);
            fratesem = FRsemAll(inds);
            
            if max(rendsinbout)<minNumRends
                continue
            end
            
            if isempty(ffmeans);
                continue
            end
            % subtract mean FF of first rend
            ffmeans = ffmeans-ffmeans(1);
            
            % === PLOT
            shadedErrorBar(rendsinbout, ffmeans, ffsems, {'Color', plotcols{j}}, 1);
        end
        lt_plot_zeroline;
        
        % FIRING
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i) ' neuron' num2str(ii)]);
        
        plotcols = lt_make_plot_colors(max(MotifNumAll), 0, 0);
        for j=1:max(MotifNumAll)
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            disp([num2str(RendInBoutAll(inds)) '-' num2str(MotifNumAll(inds))]);
            
            rendsinbout = RendInBoutAll(inds);
            motifnum = MotifNumAll(inds);
            ffmeans = FFmeanAll(inds);
            ffsems = FFsemAll(inds);
            fratemeans = FRmeanAll(inds);
            fratesem = FRsemAll(inds);
            
            if max(rendsinbout)<minNumRends
                continue
            end
            
            if isempty(ffmeans);
                continue
            end
            % subtract mean FF of first rend
            fratemeans = fratemeans-fratemeans(1);
            
            % === PLOT
            shadedErrorBar(rendsinbout, fratemeans, fratesem, {'Color', plotcols{j}}, 1);
        end
        lt_plot_zeroline;
        
        % CORRELATION
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i) ' neuron' num2str(ii)]);
        xlabel('frate (hz)'); ylabel('ff');
        plotcols = lt_make_plot_colors(max(MotifNumAll), 0, 0);
        for j=1:max(MotifNumAll)
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            disp([num2str(RendInBoutAll(inds)) '-' num2str(MotifNumAll(inds))]);
            
            rendsinbout = RendInBoutAll(inds);
            motifnum = MotifNumAll(inds);
            ffmeans = FFmeanAll(inds);
            ffsems = FFsemAll(inds);
            fratemeans = FRmeanAll(inds);
            fratesem = FRsemAll(inds);
            
            if max(rendsinbout)<minNumRends
                continue
            end
            
            if isempty(ffmeans);
                continue
            end
            fratemeans = fratemeans - fratemeans(1);
            ffmeans = ffmeans -ffmeans(1);
            
            lt_plot(fratemeans, ffmeans, {'Color', plotcols{j}});
            
        end
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
        lt_plot_line_xy;
        
    end
end


%% ========== COLLECT ALL DATA - PLOT 1) PITCH, 2) FR, 3) CORRELATION
% IN PROGRESS!!!!
% [ALL RENDS] [QUESTION: pitch corr with LMAN activity, and how much of
% that is explained by position?]


BirdnumAll=[];
NeuronNumAll = [];
MotifNumAll = [];
PitchAll = [];
NumSpksAll = [];
RendInBoutAll = [];

NumBirds = length(SummaryStruct.birds);
for i=1:NumBirds
    
    NumNeurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:NumNeurons
        
        NumMotifs = length(MOTIFSTATS.birds(i).neurons(ii).motif);
        
        for j=1:NumMotifs;
            
            SegmentsExtract = MOTIFSTATS.birds(i).neurons(ii).motif(j).SegmentsExtract;
            
            % ==== EXTRACT STUFF
            % -- 1) ff
            if ~isfield(SegmentsExtract, 'FF_val')
                continue
            end
            FFvals = [SegmentsExtract.FF_val];
            
            if all(isnan(FFvals))
                continue
            end
            
            % -- 1b) rend pos
            RendInBout = [SegmentsExtract.BOUT_RendInBout];
            
            % 2) FR
            numtrials = length(SegmentsExtract);
            clustnum = MOTIFSTATS.birds(i).neurons(ii).clustnum;
            NumSpks = [];
            for k=1:numtrials
                inds=SegmentsExtract(k).spk_Clust==clustnum;
                spktimes=SegmentsExtract(k).(spktimefield)(inds);
                
                % look at a small time window
                windowtmp_forspks = [motif_predur+window_meanFF(1) motif_predur+window_meanFF(2)];
                spktimes = spktimes(spktimes>windowtmp_forspks(1) & spktimes<windowtmp_forspks(2));
                
                NumSpks = [NumSpks length(spktimes)];
            end
            
            % =============== OUTPUT
            N = length(FFvals);
            
            BirdnumAll=[BirdnumAll i*ones(1, N)];
            NeuronNumAll = [NeuronNumAll ii*ones(1,N)];
            MotifNumAll = [MotifNumAll j*ones(1,N)];
            
            PitchAll = [PitchAll FFvals];
            NumSpksAll = [NumSpksAll NumSpks];
            RendInBoutAll = [RendInBoutAll RendInBout];
            
        end
    end
end


%% ============= PLOT [REGRESSION, PREDICTING PITCH]

count = 1;
lt_figure; hold on;

XTickLab = {};
for i=1:max(BirdnumAll)
    
    for ii=1:max(NeuronNumAll)
        
        plotcollist = lt_make_plot_colors(max(MotifNumAll), 0,0);
        % --- PITCH
        for j=1:max(MotifNumAll)
            plotcol = plotcollist{j};
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            
            rendsinbout = RendInBoutAll(inds);
            ffvals = PitchAll(inds);
            numspks = NumSpksAll(inds);
            
            
            if max(rendsinbout)<2
                continue
            end
            
            if isempty(ffvals);
                continue
            end
            
            
            % ====== REGRESSION [MULTIPLE REGR]
            % -- predict ffvals based on 1) rendition and 2) numspks
            N = length(ffvals);
            X = [ones(N,1) rendsinbout' numspks'];
            y = ffvals';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            % === PLOT
            % 1) each b by bird
            lt_subplot(4,1,1); hold on; % effect of rends
            title('effect of position (mult)');
            ylabel('(beta, predicting pitch)')
            if bint(2,1)*bint(2,2)>0 % then significant
                lt_plot(count, b(2), {'Color', plotcol});
            else
                plot(count, b(2), 'o','Color', plotcol);
            end
            line([count count], bint(2,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            lt_subplot(4,1,2); hold on; % effect of spks
            title('effect of spks (mult)');
            if bint(3,1)*bint(3,2)>0
                lt_plot(count, b(3), {'Color', plotcol});
            else
                plot(count, b(3), 'o', 'Color', plotcol);
            end
            line([count count], bint(3,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            
            % ========= REGRESSION [posotion, SINGLE]
            % 1) --- position
            % -- predict ffvals based on 1) rendition and 2) numspks
            N = length(ffvals);
            X = [ones(N,1) rendsinbout'];
            y = ffvals';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            % === PLOT
            % 1) each b by bird
            lt_subplot(4,1,3); hold on; % effect of rends
            title('effect of position');
            if bint(2,1)*bint(2,2)>0 % then significant
                lt_plot(count, b(2), {'Color', plotcol});
            else
                plot(count, b(2), 'o','Color', plotcol);
            end
            line([count count], bint(2,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            % ========= REGRESSION [frate, SINGLE]
            % 1) --- position
            % -- predict ffvals based on 1) rendition and 2) numspks
            N = length(ffvals);
            X = [ones(N,1) numspks'];
            y = ffvals';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            % === PLOT
            % 1) each b by bird
            lt_subplot(4,1,4); hold on; % effect of rends
            title('effect of spks');
            if bint(2,1)*bint(2,2)>0 % then significant
                lt_plot(count, b(2), {'Color', plotcol});
            else
                plot(count, b(2), 'o','Color', plotcol);
            end
            line([count count], bint(2,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            
            %
            % -------------------------
            count = count+1;
            XTickLab = [XTickLab [num2str(i) '-' num2str(ii)]];
        end
    end
end
lt_subplot(4,1,4);
set(gca, 'XTick', 1:length(XTickLab));
set(gca, 'XTickLabel', XTickLab);
rotateXLabels(gca, 90);
xlabel('bird-neur');




%% ============= PLOT [REGRESSION, PREDICTING FR]

count = 1;
lt_figure; hold on;

XTickLab = {};
for i=1:max(BirdnumAll)
    
    for ii=1:max(NeuronNumAll)
        
        plotcollist = lt_make_plot_colors(max(MotifNumAll), 0,0);
        % --- PITCH
        for j=1:max(MotifNumAll)
            plotcol = plotcollist{j};
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            
            rendsinbout = RendInBoutAll(inds);
            ffvals = PitchAll(inds);
            numspks = NumSpksAll(inds);
            
            
            if max(rendsinbout)<2
                continue
            end
            
            if isempty(ffvals);
                continue
            end
            
            
            % ====== REGRESSION [MULTIPLE REGR]
            % -- predict ffvals based on 1) rendition and 2) numspks
            N = length(ffvals);
            X = [ones(N,1) rendsinbout' ffvals'];
            y = numspks';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            % === PLOT
            % 1) each b by bird
            lt_subplot(4,1,1); hold on; % effect of rends
            title('effect of position (mult)');
            ylabel('(beta, predicting num spikes)')
            if bint(2,1)*bint(2,2)>0 % then significant
                lt_plot(count, b(2), {'Color', plotcol});
            else
                plot(count, b(2), 'o','Color', plotcol);
            end
            line([count count], bint(2,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            lt_subplot(4,1,2); hold on; % effect of spks
            title('effect of pitch (mult)');
            if bint(3,1)*bint(3,2)>0
                lt_plot(count, b(3), {'Color', plotcol});
            else
                plot(count, b(3), 'o', 'Color', plotcol);
            end
            line([count count], bint(3,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            
            % ========= REGRESSION [posotion, SINGLE]
            % 1) --- position
            % -- predict ffvals based on 1) rendition and 2) numspks
            N = length(ffvals);
            X = [ones(N,1) rendsinbout'];
            y = numspks';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            % === PLOT
            % 1) each b by bird
            lt_subplot(4,1,3); hold on; % effect of rends
            title('effect of position');
            if bint(2,1)*bint(2,2)>0 % then significant
                lt_plot(count, b(2), {'Color', plotcol});
            else
                plot(count, b(2), 'o','Color', plotcol);
            end
            line([count count], bint(2,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            % ========= REGRESSION [frate, SINGLE]
            % 1) --- position
            % -- predict ffvals based on 1) rendition and 2) numspks
            N = length(ffvals);
            X = [ones(N,1) ffvals'];
            y = numspks';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            % === PLOT
            % 1) each b by bird
            lt_subplot(4,1,4); hold on; % effect of rends
            title('effect of pitch');
            if bint(2,1)*bint(2,2)>0 % then significant
                lt_plot(count, b(2), {'Color', plotcol});
            else
                plot(count, b(2), 'o','Color', plotcol);
            end
            line([count count], bint(2,:), 'Color', plotcol);
            lt_plot_zeroline;
            
            
            %
            % -------------------------
            count = count+1;
            XTickLab = [XTickLab [num2str(i) '-' num2str(ii)]];
        end
    end
end
lt_subplot(4,1,4);
set(gca, 'XTick', 1:length(XTickLab));
set(gca, 'XTickLabel', XTickLab);
rotateXLabels(gca, 90);
xlabel('bird-neur');

%% ====== HISTOGRAMS (EFFECT OF POS on PITCH/FRATE)

Beta_Spikes = [];
Beta_Spikes_CI = [];
Beta_Spikes_significant = [];

Beta_Pitch = [];
Beta_Pitch_CI = [];
Beta_Pitch_significant = [];

for i=1:max(BirdnumAll)
    
    for ii=1:max(NeuronNumAll)
        
        % --- PITCH
        for j=1:max(MotifNumAll)
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            
            rendsinbout = RendInBoutAll(inds);
            ffvals = PitchAll(inds);
            numspks = NumSpksAll(inds);
            
            
            if max(rendsinbout)<2
                continue
            end
            
            if isempty(ffvals);
                continue
            end
            
            
            % ========= REGRESSION (spikes ~ position)
            % 1) --- SPIKES
            N = length(ffvals);
            X = [ones(N,1) rendsinbout'];
            y = numspks';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            Beta_Spikes = [Beta_Spikes b(2)];
            Beta_Spikes_CI = [Beta_Spikes_CI; bint(2,:)];
            if bint(2,1)*bint(2,2)>0
                Beta_Spikes_significant = [Beta_Spikes_significant 1];
            else
                Beta_Spikes_significant = [Beta_Spikes_significant 0];
            end
            
            % ---- 2) PiTCH
            N = length(ffvals);
            X = [ones(N,1) rendsinbout'];
            y = ffvals';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            Beta_Pitch = [Beta_Pitch b(2)];
            Beta_Pitch_CI = [Beta_Pitch_CI; bint(2,:)];
            if bint(2,1)*bint(2,2)>0
                Beta_Pitch_significant = [Beta_Pitch_significant 1];
            else
                Beta_Pitch_significant = [Beta_Pitch_significant 0];
            end
            
        end
    end
end

% ========= PLOT
lt_figure; hold on;

% lt_subplot(3,1,1); hold on;
xlabel('beta (spikes ~ position)', 'Color', 'b');
ylabel('beta (pitch ~ position)', 'Color', 'm');

for i=1:length(Beta_Spikes)
    if Beta_Pitch_significant(i)==1 & Beta_Spikes_significant(i)==1
        plot(Beta_Spikes(i), Beta_Pitch(i), 'ko', 'MarkerFaceColor', 'k');
    elseif Beta_Pitch_significant(i)==1 & Beta_Spikes_significant(i)==0
        plot(Beta_Spikes(i), Beta_Pitch(i), 'ko', 'MarkerFaceColor', 'm');
    elseif Beta_Pitch_significant(i)==0 & Beta_Spikes_significant(i)==1
        plot(Beta_Spikes(i), Beta_Pitch(i), 'ko', 'MarkerFaceColor', 'b');
    else
        plot(Beta_Spikes(i), Beta_Pitch(i), 'ko');
    end
    
    % -- plot CI
    %     line(Beta_Spikes_CI(i, :), [Beta_Pitch(i) Beta_Pitch(i)], 'Color','m');
    %     line([Beta_Spikes(i) Beta_Spikes(i)], Beta_Pitch_CI(i, :), 'Color','m');
end
lt_plot_zeroline;
lt_plot_zeroline_vert;

% %  =--2) historgram
% lt_subplot(3,1,2); hold on;
% title('beta (spikes ~ position)')
% lt_plot_


%% ====== HISTOGRAMS (EFFECT OF POS on PITCH/FRATE)
% FOR PITCH VS. POSITION, COMBINING ACROSS ALL NEURONS (making sure to
% avoid pseudoreplication)

Beta_Spikes = [];
Beta_Spikes_CI = [];
Beta_Spikes_significant = [];

Beta_Pitch = [];
Beta_Pitch_CI = [];
Beta_Pitch_significant = [];

for i=1:max(BirdnumAll)
    
    %% for each motif, combine all pitch data to get corr with position
% CONFIRMED BY hand thru some examples that this works
    BetaPitch_NeuronsCombined = nan(max(MotifNumAll), 1);
    BetaPitchCI_NeuronsCombined = nan(max(MotifNumAll), 2);
    
    for j=1:max(MotifNumAll)
        inds = BirdnumAll==i & MotifNumAll==j;
        ffvals = PitchAll(inds);
        rendsinbout = RendInBoutAll(inds);
        
        if max(rendsinbout)<2
            continue
        end
        if isempty(ffvals)
            continue
        end
        
        [~, inds2] = unique(ffvals); % make sure don't pseudoreplicate if got simultaneous neurons
        if length(inds2)<length(ffvals)
            disp(['bird ' num2str(i) ', motif ' num2str(j) ': kept n=' num2str(length(inds2)) ' from ' num2str(length(ffvals))]);
        else
            disp(['bird ' num2str(i) ', motif ' num2str(j) ': kept all']);
        end
        ffvals = ffvals(inds2);
        rendsinbout = rendsinbout(inds2);
        
        
        % === regression
        N = length(ffvals);
        X = [ones(N,1) rendsinbout'];
        y = ffvals';
        
        [b, bint, ~, ~, stats] = regress(y, X);
        
        BetaPitch_NeuronsCombined(j) = b(2);
        BetaPitchCI_NeuronsCombined(j,:) = bint(2,:);
        
    end
    
    %%
    for ii=1:max(NeuronNumAll)
        
        % --- PITCH
        for j=1:max(MotifNumAll)
            
            inds = BirdnumAll==i & NeuronNumAll==ii & MotifNumAll==j;
            
            rendsinbout = RendInBoutAll(inds);
            ffvals = PitchAll(inds);
            numspks = NumSpksAll(inds);
            
            
            if max(rendsinbout)<2
                continue
            end
            
            if isempty(ffvals);
                continue
            end
            
            
            % ========= REGRESSION (spikes ~ position)
            % 1) --- SPIKES
            N = length(ffvals);
            X = [ones(N,1) rendsinbout'];
            y = numspks';
            
            [b, bint, ~, ~, stats] = regress(y, X);
            
            Beta_Spikes = [Beta_Spikes b(2)];
            Beta_Spikes_CI = [Beta_Spikes_CI; bint(2,:)];
            if bint(2,1)*bint(2,2)>0
                Beta_Spikes_significant = [Beta_Spikes_significant 1];
            else
                Beta_Spikes_significant = [Beta_Spikes_significant 0];
            end
            
            % ---- 2) PiTCH
            % use the value from data combined across neurons
            Beta_Pitch = [Beta_Pitch BetaPitch_NeuronsCombined(j)];
            Beta_Pitch_CI = [Beta_Pitch_CI BetaPitchCI_NeuronsCombined(j,:)];
            
            if BetaPitchCI_NeuronsCombined(j,1)*BetaPitchCI_NeuronsCombined(j,2)>0
                Beta_Pitch_significant = [Beta_Pitch_significant 1];
            else
                Beta_Pitch_significant = [Beta_Pitch_significant 0];
            end
            
        end
    end
end

% ========= PLOT
lt_figure; hold on;

% lt_subplot(3,1,1); hold on;
ylabel('beta (spikes ~ position)', 'Color', 'b');
xlabel('beta (pitch ~ position)', 'Color', 'm');

for i=1:length(Beta_Spikes)
    if Beta_Pitch_significant(i)==1 & Beta_Spikes_significant(i)==1
        %         plot(Beta_Spikes(i), Beta_Pitch(i), 'ko', 'MarkerFaceColor', 'k');
        plot(Beta_Pitch(i), Beta_Spikes(i), 'ko', 'MarkerFaceColor', 'g');
    elseif Beta_Pitch_significant(i)==1 & Beta_Spikes_significant(i)==0
        %         plot(Beta_Spikes(i), Beta_Pitch(i), 'ko', 'MarkerFaceColor', 'm');
        plot(Beta_Pitch(i), Beta_Spikes(i), 'ko', 'MarkerFaceColor', 'm');
    elseif Beta_Pitch_significant(i)==0 & Beta_Spikes_significant(i)==1
        %         plot(Beta_Spikes(i), Beta_Pitch(i), 'ko', 'MarkerFaceColor', 'b');
        plot(Beta_Pitch(i), Beta_Spikes(i), 'ko', 'MarkerFaceColor', 'b');
    else
        %         plot(Beta_Spikes(i), Beta_Pitch(i), 'ko');
        plot(Beta_Pitch(i), Beta_Spikes(i), 'ko');
    end
    
    
    % -- plot CI
    %     line(Beta_Spikes_CI(i, :), [Beta_Pitch(i) Beta_Pitch(i)], 'Color','m');
    %     line([Beta_Spikes(i) Beta_Spikes(i)], Beta_Pitch_CI(i, :), 'Color','m');
end
% regression
lt_regress(Beta_Spikes, Beta_Pitch, 0, 0, 1, 1, 'r', 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% %  =--2) historgram
% lt_subplot(3,1,2); hold on;
% title('beta (spikes ~ position)')
% lt_plot_


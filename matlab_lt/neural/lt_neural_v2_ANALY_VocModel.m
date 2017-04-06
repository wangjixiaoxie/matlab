function VOCALSTRUCTall = lt_neural_v2_ANALY_VocModel(SummaryStruct, Binparams, plotRaw, ConvOrDiv, saveOn)
%% LT, 3/22/17 - use features for each vocalization to predict FR at various time bins

% plotRaw = 0; % individaul neuronf igs (note, hand entering neuron ind currently, in code)
% Binparams.Pretime = 0.5; % to start getting binned data (rel to onset)
% Binparams.Posttime = 0.5; % to stop getting dat (rel to onset);
% Binparams.Binsize = 0.025; % for getting spike counts

do2by2 = 0; % then uses data for bk7, where is a clear 2 x 2 and thus can do anova interaction.

%%

Pretime = Binparams.Pretime;
Posttime = Binparams.Posttime;
Binsize = Binparams.Binsize;

binedges = 0:Binsize:(Pretime+Posttime);
mindattoplot = 5; % min N for a given category.

assert(mod(Pretime, Binsize)==0 && mod(Posttime, Binsize)==0, 'pre/post time not mulitple of binsize...');


%%

VOCALSTRUCTall = struct;
VOCALSTRUCTall.global.binedges = binedges;
VOCALSTRUCTall.global.pretime = Pretime;
VOCALSTRUCTall.global.posttime = Posttime;
VOCALSTRUCTall.global.bincenters = binedges+Binsize/2;
VOCALSTRUCTall.global.mindattoplot = mindattoplot;

Numbirds = length(SummaryStruct.birds);
for z = 1:Numbirds
    Numneurons = length(SummaryStruct.birds(z).neurons);
    for zz = 1:Numneurons
        
        datstruct = SummaryStruct.birds(z).neurons(zz);
        
        % ==================== EXTRACT
        cd(datstruct.dirname);
        
        % load dat for neuron
        cd ..
        batchf = datstruct.batchfilename;
        chan = datstruct.channel;
        [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, chan);
        
        %%
        
        % go thru all vocalizations, for each, extract data
        tic
        counter = 0;
        numsyls = length(SongDat.AllLabels);
        VOCALSTRUCT = struct;
        for j=1:numsyls
            syl=SongDat.AllLabels(j);
            
            % ------ skip if is blank
            if strcmp(syl, '-')
                continue
            end
            
            % ----- sequence information
            if j>10
                pre_ten_syls=SongDat.AllLabels(j-10:j-1);
            elseif j>1
                pre_ten_syls=SongDat.AllLabels(1:j-1);
            else
                pre_ten_syls = '?';
            end
            if length(SongDat.AllLabels)>j+9
                post_ten_syls=SongDat.AllLabels(j+1:j+10);
            elseif length(SongDat.AllLabels)>j
                post_ten_syls=SongDat.AllLabels(j+1:end);
            else
                post_ten_syls = '?';
            end
            
            presyl = pre_ten_syls(end);
            postsyl = post_ten_syls(1);
            
            % -- syl dur
            Syldur = SongDat.AllOffsets(j)-SongDat.AllOnsets(j);
            
            % ---- gap information
            Gapdur_post = [];
            Gapdur_pre = [];
            try
                Gapdur_post = SongDat.AllOnsets(j+1) - SongDat.AllOffsets(j);
            catch err
            end
            try
                Gapdur_pre = SongDat.AllOnsets(j) - SongDat.AllOffsets(j-1);
            catch err
            end
            
            % --- song file/time information
            fs=NeurDat.metaDat(1).fs;
            globalOnsetTime=SongDat.AllOnsets(j); % sec
            globalOffsetTime = SongDat.AllOffsets(j);
            
            cumOnsSamps=cumsum([0 NeurDat.metaDat.numSamps]);
            globalOnsetSamp=globalOnsetTime*fs;
            songind=find((globalOnsetSamp-cumOnsSamps)>0, 1, 'last');
            songfname=NeurDat.metaDat(songind).filename;
            
            % ----- spiking information, multiple time bins
            if strcmp(ConvOrDiv, 'conv')
                onsettmp = 1000*(globalOnsetTime - Pretime);
                offsettmp = 1000*(globalOnsetTime + Posttime);
            elseif strcmp(ConvOrDiv, 'div');
                onsettmp = 1000*(globalOffsetTime - Pretime);
                offsettmp = 1000*(globalOffsetTime + Posttime);
            end
            
            inds = NeurDat.spikes_cat.cluster_class(:,2) > onsettmp ...
                & NeurDat.spikes_cat.cluster_class(:,2) < offsettmp ...
                & NeurDat.spikes_cat.cluster_class(:,1) == datstruct.clustnum;
            Spktimes_relonset = NeurDat.spikes_cat.cluster_class(inds, 2) - ...
                onsettmp;
            Spktimes_relonset = Spktimes_relonset/1000; % seconds
            
            % bin the spike times
            [Bincounts] = histc(Spktimes_relonset, binedges);
            % --- troubleshooting using old method
            % Bincounts = zeros(1, 13);
            % Yspks = {Spktimes_relonset};
            % window_fr = 0.025;
            % windshift_fr = 0.025;
            %         [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
            %         lt_neural_plotRastMean(Yspks, window_fr, windshift_fr, 0, '');
            % Bincounts(1:length(ymean_hz)) = ymean_hz;
            
            
            % =================== OUTPUT
            counter = counter +1;
            VOCALSTRUCT.data_vocalrend(counter).syl = syl;
            VOCALSTRUCT.data_vocalrend(counter).presyl = presyl;
            VOCALSTRUCT.data_vocalrend(counter).postsyl = postsyl;
            
            VOCALSTRUCT.data_vocalrend(counter).pre_ten_syls = pre_ten_syls;
            VOCALSTRUCT.data_vocalrend(counter).post_ten_syls = post_ten_syls;
            VOCALSTRUCT.data_vocalrend(counter).syldur = Syldur;
            VOCALSTRUCT.data_vocalrend(counter).Gapdur_post = Gapdur_post;
            VOCALSTRUCT.data_vocalrend(counter).Gapdur_pre = Gapdur_pre;
            VOCALSTRUCT.data_vocalrend(counter).globalOnsetTime = globalOnsetTime;
            VOCALSTRUCT.data_vocalrend(counter).globalOffsetTime = globalOffsetTime;
            VOCALSTRUCT.data_vocalrend(counter).songind = songind;
            VOCALSTRUCT.data_vocalrend(counter).songfname = songfname;
            VOCALSTRUCT.data_vocalrend(counter).Spktimes_relonset = Spktimes_relonset;
            if size(Bincounts,2)>1
                Bincounts= Bincounts';
            end
            VOCALSTRUCT.data_vocalrend(counter).Spk_bincounts = Bincounts';
            
            
            
            VOCALSTRUCT.global.fs = fs;
        end
        
        
        % ===== COMBINE ACROSS NERUONS/BIRDS
        VOCALSTRUCTall.birds(z).neurons(zz).dat = VOCALSTRUCT;
        
        toc
    end
end


% if (0)
%     %% ========================= TEMPORARY, FILTER VOCALSTRUCT
%     % 2 x 2 transition matrix (full) for bk7
%     % IMPORTANT: must change ii for loop to presyltypes during figures;
%     
%     if do2by2 ==1
%         
%         PresylPool = {'p','h'};
%         SylPool = {'p','r'};
%         DesiredBird = 'bk7';
%         useinteraction=1;
%         
%         for z = 1:Numbirds
%             
%             if ~strcmp(SummaryStruct.birds(z).birdname, DesiredBird);
%                 continue
%             end
%             
%             Numneurons = length(SummaryStruct.birds(z).neurons);
%             for zz = 1:Numneurons
%                 
%                 Vocalstruct = VOCALSTRUCTall.birds(z).neurons(zz).dat;
%                 
%                 inds = [];
%                 for j=1:length(Vocalstruct.data_vocalrend)
%                     
%                     if any(strcmp(Vocalstruct.data_vocalrend(j).syl, SylPool)) ...
%                             & any(strcmp(Vocalstruct.data_vocalrend(j).presyl, PresylPool));
%                         inds = [inds j];
%                     end
%                 end
%                 
%                 % resave output
%                 VOCALSTRUCTall.birds(z).neurons(zz).dat.data_vocalrend = Vocalstruct.data_vocalrend(inds);
%             end
%         end
%     else
%         useinteraction=0;
%     end
%     
%     
%     
%     %% ==== PLOT
%     z=1; % pick out one to plot
%     zz=12;
%     timebin = 27;
%     if plotRaw==1
%         
%         VOCALSTRUCT = VOCALSTRUCTall.birds(z).neurons(zz).dat;
%         
%         lt_figure; hold on;
%         plotMean = 0;
%         
%         % ==== 1) for single time bin, single neuron, plot data
%         lt_subplot(2,1,1); hold on;
%         tmp1 = binedges(timebin);
%         tmp2 = binedges(timebin+1);
%         title(['ignoring context; timebin: ' num2str(tmp1) ' to ' num2str(tmp2) ' (onset=' num2str(Pretime) ') sec']);
%         ylabel('spks/bin');
%         
%         % x axis for syl, vary color for sequence (1-back)
%         syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
%         for i=1:length(syltypes)
%             syl = syltypes{i};
%             inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
%             
%             if sum(inds)<mindattoplot
%                 continue % not enough samples
%             end
%             
%             tmp = {VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts};
%             SpkCntMat = cell2mat(tmp');
%             
%             % -- extract each trial dat
%             Y = SpkCntMat(:, timebin);
%             X = i;
%             
%             if plotMean ==1
%                 lt_plot(X, mean(Y), {'Errors', lt_sem(Y), 'Color', 'k'});
%             else
%                 distributionPlot(Y, 'xValues', X, 'addSpread', 0, 'showMM', 0);
%             end
%             %     plot(X, Y, 'o');
%         end
%         lt_plot_zeroline;
%         
%         
%         
%         % ==== 2) for single time bin, single neuron, plot data [including sequence
%         % info]
%         lt_subplot(2,1,2); hold on;
%         title('colors = context');
%         xlabel('syl');
%         ylabel('spks/bin');
%         grid on
%         % x axis for syl, vary color for sequence (1-back)
%         syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
%         plotcols = lt_make_plot_colors(length(syltypes), 0, 0);
%         
%         for i=1:length(syltypes)
%             syl = syltypes{i};
%             %     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
%             
%             % --- break out into preceding syl
%             %         presyltypes = unique({VOCALSTRUCT.data_vocalrend.presyl});
%             for ii=1:length(syltypes)
%                 sylpre = syltypes{ii};
%                 
%                 inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl) & ...
%                     strcmp({VOCALSTRUCT.data_vocalrend.presyl}, sylpre);
%                 
%                 if sum(inds)<mindattoplot
%                     disp('NO')
%                     continue % not enough samples
%                 else
%                     disp('YES')
%                 end
%                 
%                 tmp = {VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts};
%                 SpkCntMat = cell2mat(tmp');
%                 
%                 % -- extract each trial dat
%                 Y = SpkCntMat(:, timebin);
%                 X = i + 0.3-0.6*rand;
%                 
%                 if plotMean==1
%                     lt_plot(X, mean(Y), {'Errors', lt_sem(Y), 'Color', plotcols{ii}});
%                     lt_plot_text(X+0.01, mean(Y), sylpre, plotcols{ii});
%                 else
%                     distributionPlot(Y, 'xValues', X, 'addSpread', 0, 'color', plotcols{ii}, ...
%                         'distWidth', 0.3, 'showMM', 0);
%                     %     distributionPlot(Y, 'xValues', X, 'addSpread', 0, 'color', plotcols{ii});
%                 end
%                 %     plot(X, Y, 'o');
%             end
%         end
%         lt_plot_zeroline;
%         set(gca, 'XTick', 1:length(syltypes), 'XTickLabel', syltypes);
%         
%         %% DIAGNOSTIC (WAS COMPARING IT TO FORMER CODE, M AKING SURE WAS SIMILAR
%         % =============== PLOT RUNNING TIMECOURSE [no smoothuing]
%         figcount=1;
%         subplotrows=8;
%         subplotcols=3;
%         fignums_alreadyused=[];
%         hfigs=[];
%         
%         % x axis for syl, vary color for sequence (1-back)
%         syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
%         plotcols = lt_make_plot_colors(length(syltypes), 0, 0);
%         
%         for i=1:length(syltypes)
%             syl = syltypes{i};
%             %     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
%             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%             title(['to ' syl]);
%             
%             % --- break out into preceding syl
%             presyltypes = unique({VOCALSTRUCT.data_vocalrend.presyl});
%             for ii=1:length(syltypes)
%                 sylpre = syltypes{ii};
%                 
%                 inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl) & ...
%                     strcmp({VOCALSTRUCT.data_vocalrend.presyl}, sylpre);
%                 
%                 if sum(inds)<mindattoplot
%                     continue % not enough samples
%                 end
%                 
%                 if (0)
%                     SpkCntMat = reshape([VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts], sum(inds), ...
%                         length(VOCALSTRUCTall.global.binedges)); % collect all into matrix of spk counts (N x bins)
%                 else
%                     tmp = {VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts};
%                     SpkCntMat = cell2mat(tmp');
%                 end
%                 
%                 % -- extract each trial dat
%                 Y = mean(SpkCntMat, 1);
%                 X = VOCALSTRUCTall.global.binedges(1:end)+Binsize/2;
%                 
%                 plot(X, Y./Binsize, 'Color', plotcols{ii});
%                 
%                 line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim ,'Color', 'k')
%             end
%         end
%         lt_plot_zeroline;
%         set(gca, 'XTick', 1:length(syltypes), 'XTickLabel', syltypes);
%         
%         
%         % =============== PLOT RUNNING TIMECOURSE [smoothing]
%         if (0)
%             figcount=1;
%             subplotrows=8;
%             subplotcols=3;
%             fignums_alreadyused=[];
%             hfigs=[];
%             
%             % x axis for syl, vary color for sequence (1-back)
%             syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
%             plotcols = lt_make_plot_colors(length(syltypes), 0, 0);
%             
%             for i=1:length(syltypes)
%                 syl = syltypes{i};
%                 %     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
%                 [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%                 title(['to ' syl]);
%                 
%                 % --- break out into preceding syl
%                 %             presyltypes = unique({VOCALSTRUCT.data_vocalrend.presyl});
%                 for ii=1:length(syltypes)
%                     sylpre = syltypes{ii};
%                     
%                     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl) & ...
%                         strcmp({VOCALSTRUCT.data_vocalrend.presyl}, sylpre);
%                     
%                     if sum(inds)<mindattoplot
%                         continue % not enough samples
%                     end
%                     
%                     Yspks = {VOCALSTRUCT.data_vocalrend(inds).Spktimes_relonset};
%                     
%                     window_fr = 0.04;
%                     windshift_fr = 0.002;
%                     [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
%                         lt_neural_plotRastMean(Yspks, window_fr, windshift_fr, 0, '');
%                     
%                     shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{ii}}, 1);
%                     
%                     line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim ,'Color', 'k')
%                 end
%             end
%             lt_plot_zeroline;
%             set(gca, 'XTick', 1:length(syltypes), 'XTickLabel', syltypes);
%         end
%     end
%     
%     
%     
%     %% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ ANOVA
%     %% ==================== EXTRACT DATA
%     % useinteraction =1; use for if do 2 x 2 above.
%     if do2by2==1
%         Groupnames = {'Syl', 'Presyl'}; % adds interaction, and pulls out stats assuming there's only 2 factors
%     else
%         Groupnames = {'Syl', 'Presyl', 'Postsyl'};
%     end
%     numtimebins = length(VOCALSTRUCTall.global.bincenters);
%     ANOVAOUTPUT = struct;
%     logtransform = 1;
%     
%     % lt_figure; hold on;
%     % title('anova rel to syl onset');
%     Numbirds = length(VOCALSTRUCTall.birds);
%     
%     for z=1:Numbirds
%         Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
%         %     z=1
%         for zz=1:Numneurons
%             VOCALSTRUCT = VOCALSTRUCTall.birds(z).neurons(zz).dat;
%             %         zz=12
%             for i=1:numtimebins
%                 
%                 % ==========
%                 tmp = {VOCALSTRUCT.data_vocalrend.Spk_bincounts};
%                 SpkCntMat = cell2mat(tmp');
%                 SpkCnts = SpkCntMat(:, i);
%                 
%                 if logtransform==1
%                     SpkCnts = log10(SpkCnts+1);
%                 end
%                 
%                 Syl = [VOCALSTRUCT.data_vocalrend.syl]';
%                 Presyl = [VOCALSTRUCT.data_vocalrend.presyl]';
%                 Postsyl = [VOCALSTRUCT.data_vocalrend.postsyl]';
%                 
%                 % ------ 1) LME model
%                 if (0)
%                     modelform = 'SpkCnts ~ 1 + Syl + Presyl';
%                     X = table(SpkCnts, Syl, Presyl);
%                     lme = fitlme(X, modelform);
%                 end
%                 
%                 % ----- 2) ANOVA
%                 Group = {};
%                 for j=1:length(Groupnames)
%                     Group{j} = eval(Groupnames{j});
%                 end
%                 if useinteraction ==1;
%                     [p, tabletmp] = anovan(SpkCnts, Group, 'display', 'off', ...
%                         'model', 'interaction');
%                 else
%                     [p, tabletmp] = anovan(SpkCnts, Group, 'display', 'off');
%                 end
%                 
%                 numgroups = length(Group)+useinteraction;
%                 ss_total = tabletmp{3+numgroups,2};
%                 ms_error = tabletmp{2+numgroups,5};
%                 
%                 Eta2_all = [];
%                 Omega2_all = [];
%                 
%                 for j=1:numgroups
%                     ss_effect = tabletmp{1+j,2};
%                     df_effect = tabletmp{1+j,3};
%                     
%                     numerator = ss_effect - df_effect*ms_error;
%                     denominator = ss_total + ms_error;
%                     omega2 = numerator/denominator;
%                     eta2 = ss_effect/ss_total;
%                     
%                     Eta2_all = [Eta2_all eta2];
%                     Omega2_all = [Omega2_all omega2];
%                 end
%                 
%                 ANOVAOUTPUT.timebin(i).Eta2_all = Eta2_all;
%                 ANOVAOUTPUT.timebin(i).Omega2_all = Omega2_all;
%             end
%             ANOVAOUTPUT.global.Groupnames = Groupnames;
%             if useinteraction ==1
%                 ANOVAOUTPUT.global.Groupnames = [Groupnames 'interaction'];
%             end
%             VOCALSTRUCTall.birds(z).neurons(zz).anova.output = ANOVAOUTPUT;
%             
%             if (0)
%                 % =================== PLOT ANOVA RESULTS
%                 tmp = {ANOVAOUTPUT.timebin.Omega2_all};
%                 OmegaAcrossBins = cell2mat(tmp');
%                 
%                 numgroups = length(ANOVAOUTPUT.global.Groupnames);
%                 plotcols = lt_make_plot_colors(length(Groupnames), 0,0);
%                 for i=1:numgroups
%                     groupname = ANOVAOUTPUT.global.Groupnames{i};
%                     
%                     omegavals = OmegaAcrossBins(:,i);
%                     xvals = VOCALSTRUCTall.global.bincenters;
%                     
%                     plot(xvals, omegavals, '-', 'Color' ,plotcols{i});
%                 end
%             end
%         end
%     end
%     
%     
%     
%     %% ================== PLOT ANOVA RESULTS
%     figcount=1;
%     subplotrows=3;
%     subplotcols=2;
%     fignums_alreadyused=[];
%     hfigs=[];
%     
%     numgroups = length(ANOVAOUTPUT.global.Groupnames);
%     plotcols = lt_make_plot_colors(length(ANOVAOUTPUT.global.Groupnames), 0,0);
%     
%     
%     % ============ FOR A SINGLE NEURON
%     if (0)
%         z=1;
%         zz=12;
%         
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         birdname = SummaryStruct.birds(z).birdname;
%         xlabel('timebin');
%         ylabel('omega^2');
%         
%         ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
%         for g =1:numgroups
%             groupname = Groupnames{g};
%             
%             % =================== PLOT ANOVA RESULTS
%             tmp = {ANOVAOUTPUT.timebin.Omega2_all};
%             OmegaAcrossBins = cell2mat(tmp');
%             
%             omegavals = OmegaAcrossBins(:,g);
%             xvals = VOCALSTRUCTall.global.bincenters;
%             
%             plot(xvals, omegavals, '-', 'Color' ,plotcols{g});
%         end
%         lt_plot_zeroline;
%         line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
%         ylim([-0.05 0.2]);
%     end
%     
%     
%     % ======== ONE FOR EACH BIRD
%     for z=1:Numbirds
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         birdname = SummaryStruct.birds(z).birdname;
%         title(['bird: ' birdname]);
%         xlabel('timebin');
%         ylabel('omega^2');
%         
%         Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
%         for zz=1:Numneurons
%             ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
%             for g =1:numgroups
%                 groupname = ANOVAOUTPUT.global.Groupnames{g};
%                 
%                 % =================== PLOT ANOVA RESULTS
%                 tmp = {ANOVAOUTPUT.timebin.Omega2_all};
%                 OmegaAcrossBins = cell2mat(tmp');
%                 
%                 omegavals = OmegaAcrossBins(:,g);
%                 xvals = VOCALSTRUCTall.global.bincenters;
%                 
%                 plot(xvals, omegavals, '-', 'Color' ,plotcols{g});
%             end
%         end
%         lt_plot_zeroline;
%         line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
%         ylim([-0.05 0.2]);
%     end
%     
%     
%     % ==== 1) ONE FOR EACH GROUP, EACH SYL ONE LINE
%     
%     for g =1 :numgroups
%         groupname = ANOVAOUTPUT.global.Groupnames{g};
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         title(groupname);
%         xlabel('timebin');
%         ylabel('omega^2');
%         
%         for z=1:Numbirds
%             Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
%             
%             for zz=1:Numneurons
%                 ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
%                 
%                 % =================== PLOT ANOVA RESULTS
%                 tmp = {ANOVAOUTPUT.timebin.Omega2_all};
%                 OmegaAcrossBins = cell2mat(tmp');
%                 
%                 omegavals = OmegaAcrossBins(:,g);
%                 xvals = VOCALSTRUCTall.global.bincenters;
%                 
%                 plot(xvals, omegavals, '-', 'Color' ,plotcols{g});
%             end
%         end
%         lt_plot_zeroline;
%         line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
%         ylim([-0.05 0.2]);
%     end
%     
%     
%     % ======= 2) combine all in one plot
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title(groupname);
%     xlabel('timebin');
%     ylabel('omega^2');
%     
%     for g =1 :numgroups
%         groupname = ANOVAOUTPUT.global.Groupnames{g};
%         
%         Omega_all = [];
%         for z=1:Numbirds
%             Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
%             
%             for zz=1:Numneurons
%                 ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
%                 
%                 % =================== PLOT ANOVA RESULTS
%                 tmp = {ANOVAOUTPUT.timebin.Omega2_all};
%                 OmegaAcrossBins = cell2mat(tmp');
%                 
%                 omegavals = OmegaAcrossBins(:,g);
%                 
%                 Omega_all = [Omega_all omegavals];
%             end
%         end
%         
%         x = VOCALSTRUCTall.global.bincenters;
%         ymean = mean(Omega_all');
%         ystd = std(Omega_all', 0,1);
%         ysem = lt_sem(Omega_all');
%         shadedErrorBar(x, ymean, ysem, {'Color', plotcols{g}},1);
%     end
%     
%     lt_plot_zeroline;
%     % legend(gca, Groupnames);
%     line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
%     axis tight
%     
%     
% end



%% ============= SAVE FIGURES
if saveOn==1
    lt_neural_v2_GoToFolder;
    
    cd FIGS/
    try
        cd('VocModel')
    catch err
        mkdir('VocModel');
        cd('VocModel');
    end
    
    tstamp = lt_get_timestamp(0);
    
    mkdir([tstamp '_' ConvOrDiv]);
    cd([tstamp '_' ConvOrDiv]);
    lt_save_all_figs
end






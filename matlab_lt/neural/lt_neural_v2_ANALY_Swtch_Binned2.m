function lt_neural_v2_ANALY_Swtch_Binned(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg, saveFigs, onlySingleDir, Bregion)
%% ------------

% binbysong = 1; % this first makes datapoints by song, not by rendition. this allows comparing targ and nontarg
% removeOutlier =1; % if 1 then removes >3std for slow change in FF and Nspks (slopes)
% useLagZeroCorr = 0; % if 0, then takes 5 trial bin centered at 0 of xcorr.
% removeTrainStart = 0; % ad hoc, remove N trials from start, i.e. exclude startle

constrainToHaveAtLeastOneSameType = 1; % if 1 then reduces min songs to make sure gets a same type 
% has to be at least 75% 

ScaleNspkByDrift = 0;
                % e.g. if drift goes down, then expect variance across syls
                % to decresae - then syls that start with lower FR would
                % show artificial increase in FR (and vice versa for those
                % starting higher) - quick solution, scale all
                % drift-subtracted values by mean FR
%%

% winsize = 29; % for regression (to get slope of learning)

assert(onlyPlotTargNontarg~=0, 'have not coded for nontarg syls yet');
assert(onlySingleDir==1, 'have not coded for when targ multiple directions yet...');

% ccmaxlagbin = 5;
% ccmaxlagtrial = 25;
% 
% if any([isempty(birdname_get) isempty(exptname_get) isempty(switchnum_get)])
%     plotraw = 0;
% else
%     plotraw =1;
% end
% 
%% PLOT FOR ALL SWITCHES


Numbirds = length(SwitchStruct.bird);
% %
% WindowToPlot = [-0.15 0.1]; % relative to syl onset, what to plot

FFsmthbinsize = 10;

% minrends = 4; % for both train and base

premotorWind = SwitchStruct.params.premotorWind;


% ------ for bins across renditions
% minbinsize = 15; % will take this or numbase rends (divbided by 2), whichever is larger
% maxbinsize = 25;
%%

% onlyPlotTargNontarg = 1; % if 1, then only targ syl, if 0, then all syls;
% if 2 then only targ


%%


if mod(FFsmthbinsize,2)==0
    FFsmthbinsize= FFsmthbinsize+1; % conver to odd, so that median tval is at an actual datapoint.
end

%%

% % -- xcorr
% AllBinned_SpkMean_vs_FFchange = [];
% AllBinned_FRsmChange_vs_FFchange = [];
% AllTrial_Spk_vs_FFslope = [];
% AllTrial_Spk_vs_FF = [];
%
% % -- corr
% AllTrial_Spk_vs_FF_RHO = [];

AllIsSame = [];
AllIsTarg = [];
AllBirdnum = [];
AllExptnum = [];
AllSwnum = [];
AllMotifnum = [];
AllNeurnum = [];

AllFFslope = [];
AllSpkSTDsongs = [];
AllFFvsSpk_corr = [];
AllFFslopeSig = [];

% AllLearnSlopeScaled = [];
% AllLearnSlopeSig = [];
%
% All_SpkMean_STDtrials = [];
% All_SpkMean_slopeScaled = [];
% %                     All_SpkMean_Change = [];
% All_SpkMean_STDtrials_notCV = [];
%
% All_SpkMean_STDtrials_base = [];
%
% All_FR_Modulation = [];
% All_StartFromWNOff = [];
% All_SingleUnit = [];

AllTrainContingency = [];

AllSpkDiff_lateminusearly = [];
AllSpkDiffFrac_lateminusearly = [];


                AllRawTvals = [];
                AllRawNspks = [];
                AllRawTrainInds = [];
        AllRawFF = [];
AllRawBaseInds = [];
AllRawNspksMinusDrift = [];
AllRawFRSmooth = {};

AllFFvsSpkminusdiff_corr = [];
AllSpkDiff_minusDiff_lateminusbase = [];
AllXcorr_NspkminusDiff_vs_FF = {};
AllPairsMotifs_FF_vs_SpkMinusDiff_corr = {};
AllBaseFFvsSpkminusdiff_corr = [];
AllBaseFFvsSpkminusdiff_corr_p = [];

for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        SummaryStruct =	MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        
        % ====== sanity check
        assert(length(SummaryStruct.birds(1).neurons) == length(MotifStats.neurons), 'asdfd');
        
        
        motiflist = MotifStats.params.motif_regexpr_str;
        targsyls = MotifStats.params.TargSyls;
        samesyls = MotifStats.params.SameTypeSyls;
        nummotifs = length(motiflist);
        
        %         WindowToPlot2 = [MotifStats.params.motif_predur+WindowToPlot(1) ...
        %             MotifStats.params.motif_predur+WindowToPlot(2)]; % rel data onset (not syl onset)
        
        for iii=1:numswitches
            
            goodneurons = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons;
            
            if isempty(goodneurons)
                disp(['---SKIPPING - ' birdname '-' exptname '-sw' num2str(iii) ' (NO GOOD NEURONS)']);
                continue
            end
            
            swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            %             plotcols = lt_make_plot_colors(max(goodneurons), 0, 0);
            
            %             if onlyPlotTargNontarg==1
            %                 motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls)) ...
            %                     find(ismember(motiflist, MotifStats.params.SameTypeSyls))];
            %             elseif onlyPlotTargNontarg==2
            %                 motifstoplot = [find(ismember(motiflist, MotifStats.params.TargSyls))];
            %             else
            %                 motifstoplot = 1:nummotifs;
            %             end
            
%             
            % --- learning at target
            targlearndir = unique(cell2mat({SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2:2:end}}));
            
            if onlySingleDir==1
                if length(targlearndir)>1
                    continue
                end
            end
            
            learnconting = unique(cell2mat(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies(2:2:end)'), 'rows');
%             if isempty(learnconting)
%                 disp('asdfdsaf');
%                 keyboard
%             end
%                 

                
%             if strcmp(birdname, 'br92br54') & strcmp(exptname, 'LMANlearn6') ...
%                     & iii==7
%                 keyboard
%             end
            
            % ==== 1) for each motif, PLOT RASTER, SMTHED, AND NEURAL/FF
            for nn=goodneurons
                
                % ============ brain region filter
                bregionThis = SummaryStruct.birds(1).neurons(nn).NOTE_Location;
                if ~isempty(Bregion)
                    if ~any(strcmp(Bregion, bregionThis))
                        continue
                    end
                end
                
                
                
                %% make matrices of (song x motif)
                figcount=1;
                subplotrows=5;
                subplotcols=2;
                fignums_alreadyused=[];
                hfigs=[];
                
                
                % ============= 1) get list of songs for each motif
                FnamesAll = cell(1, length(motiflist));
                NumSongs = nan(1, length(motiflist));
                IsTarg = nan(1, length(motiflist));
                IsSame = nan(1, length(motiflist));
                for j=1:length(motiflist)
                    
                    if ~isfield(MotifStats.neurons(nn).motif(j).SegmentsExtract, 'song_filename')
                        continue
                    end
                    
                    fnames = unique([MotifStats.neurons(nn).motif(j).SegmentsExtract.song_datenum]);
                    
                    FnamesAll{j} = fnames;
                    NumSongs(j) = length(fnames);
                    IsTarg(j) = any(strcmp(targsyls, motiflist{j}));
                    IsSame(j) = any(strcmp(samesyls, motiflist{j}));
                end
                
                if rand<0
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    lt_plot_histogram(NumSongs(~isnan(NumSongs)));
                end
                
                % ============ 2) throw out motifs that are rare
                n1 = min(1*NumSongs(IsTarg==1));
                n2 =NumSongs(IsTarg==1 | IsSame==1);
                N = min(n2(n2>=n1));

                if constrainToHaveAtLeastOneSameType==1
                n3 = NumSongs(IsSame==1); % constrain to have at least one same type (if same type exists)
                if max(n3)>0.75*N
                N = min([max(n3) N]);
                end
                end
                %                 N = min(NumSongs(IsTarg==1 | IsSame==1));
                %                 prctile(N, []
                N = round(0.95*N); % any motifs with fewer than this songs, will throw out.
                
%                 if strcmp(birdname, 'br92br54') & strcmp(exptname, 'LMANlearn4') ...
%                         & iii==3
%                    keyboard 
%                 end
%                     
                
                
                motifsToKeep= find(NumSongs>=N);
                
                commonSongs = FnamesAll{motifsToKeep(1)};
                for j=2:length(motifsToKeep)
                    commonSongs = intersect(commonSongs, FnamesAll{motifsToKeep(j)});
                end
                commonSongs = sort(commonSongs);
                
                % ============= 3) prepare matrices with songs of
                % intersection
                ThisMotifOrigNums = motifsToKeep;
                ThisNSpksMat  = nan(length(commonSongs), length(motifsToKeep));
                ThisFFmat = nan(length(commonSongs), length(motifsToKeep));
                ThisFRsmooth = cell(length(commonSongs), length(motifsToKeep));
                ThisTvalsVec = commonSongs;
                ThisTvalsVec = lt_convert_EventTimes_to_RelTimes(...
                    datestr(floor(min(ThisTvalsVec)), 'ddmmmyyyy'), ThisTvalsVec);
                ThisTvalsVec = ThisTvalsVec.FinalValue;
                ThisIsTargVec = IsTarg(motifsToKeep);
                ThisIsSameVec = IsSame(motifsToKeep);
                
                % ---- what are baseline/training songs?
                ThisbaseIndsSongs = find(commonSongs>=swpre & commonSongs<swthis);
                ThistrainIndsSongs = find(commonSongs>=swthis & commonSongs<swpost);
                
%                 
%                 if strcmp(birdname, 'wh6pk36') & strcmp(exptname, 'LMANlearn2') & ...
%                         iii==1 
%                     keyboard
%                 end
%                 
                
                
                % ============ 4) Collect data
                for j=1:length(motifsToKeep)
                    motifnum = motifsToKeep(j);
                    
                    segextract = MotifStats.neurons(nn).motif(motifnum).SegmentsExtract;
                    
                    % ============== GET SPIKES MATRIX
                    clustnum = MotifStats.neurons(nn).clustnum;
                    clustCell = {segextract.spk_Clust};
                    spkCell = {segextract.spk_Times};
                    % -- sanity check
                    assert(all(cell2mat(clustCell) == clustnum), 'asdfasd');
                    % --- extract numspikes for each trial (in premotor
                    % window)
                    windspk = MotifStats.params.motif_predur+premotorWind;
                    numtrials = length(spkCell);
                    Nspks = [];
                    for tt = 1:numtrials
                        spkt = spkCell{tt} >= windspk(1) & spkCell{tt}<windspk(2);
                        Nspks = [Nspks; length(spkt)];
                    end
                    
                    % ============= GET FFVALS MATRIX
                    % -- will subtract base and also flip if negative
                    % learnign
                    ffvals = [segextract.FF_val];
                    baseInds = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(motifnum).baseInds;
                    assert(length(baseInds) == length(ffvals), 'asdfasd');
                    ffvals = ffvals - mean(ffvals(baseInds)); % subtract baseline mean.
                    ffvals = ffvals * targlearndir; % positive = learning
                    
                    
                    % ============ GET CELL ARRAY OF SMOOTHED FR
                    % -- figure out which inds are in premotor window
                    xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    indtmp = xbin>windspk(1) & xbin<=windspk(2);
                    
                    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
                    frmat = frmat(indtmp, :);
                    
                    
                    
                    % ################## GO THRU ALL SONGS AND COMBINE DATA
                    for jj=1:length(commonSongs)
                        songthis = commonSongs(jj);
                        
                        indthis = [segextract.song_datenum] == songthis;
                        assert(any(indthis),'asdf');
                        
                        % ------------- OUTPUT
                        ThisNSpksMat(jj, j) = mean(Nspks(indthis));
                        ThisFFmat(jj, j) = mean(ffvals(indthis));
                        ThisFRsmooth{jj, j} = mean(frmat(:, indthis),2);
                        
                    end
                end
                
              %% =============== GET STATS FOR EACH MOTIF
                % -- do for each motif separately
                nmotifs = size(ThisFFmat,2);
                
                
                % ------------------ SLOPE OF LEARNING
                learnslope = nan(1,nmotifs);
                learnsig = nan(1,nmotifs);
                
                for j=1:nmotifs
                    ffvals = ThisFFmat(:,j);
                    tvals = ThisTvalsVec;
                    
                    if all(isnan(ffvals))
                        learnslope(j) = nan;
                        learnsig(j) = nan;
                        continue
                    end
                    
                    [b,bint] = ...
                        lt_regress(ffvals(ThistrainIndsSongs), ...
                        tvals(ThistrainIndsSongs), 0);
                    
                    learnslope(j) = b(2)/(bint(2,2) - bint(2,1));
                    learnsig(j) = bint(2,1)*bint(2,2)>0;
                end
                
                
                % ------------------ STD OF SPIKES
                spkstd = nan(1,nmotifs);
                for j=1:nmotifs
                    
                    nspk = ThisNSpksMat(:, j);
                    
                    spkstd(j) = std(nspk(ThistrainIndsSongs))/mean(nspk(ThistrainIndsSongs));
                end
                
                
                % ------------------ SPK vs. FF corr
                corr_FFvsSpk = nan(1,nmotifs);
                for j=1:nmotifs
                    
                    ffvals = ThisFFmat(ThistrainIndsSongs,j);
                    nspk = ThisNSpksMat(ThistrainIndsSongs,j);
                    
                    rho = corr(nspk, ffvals);
                    
                    corr_FFvsSpk(j) = rho;
                end
                
                
                % ------------------ SPK (subtract mean of diff type) vs.
                % FF corr
                indtmp = ThisIsTargVec==0 & ThisIsSameVec==0;
                assert(any(indtmp)==1, 'asfdsa');
                ThisNspksDrift = nanmean(ThisNSpksMat(:,indtmp), 2);
                ThisNSpksMat_minusDiff = ThisNSpksMat - repmat(ThisNspksDrift, 1, size(ThisNSpksMat,2));
                
                if ScaleNspkByDrift ==1
                % e.g. if drift goes down, then expect variance across syls
                % to decresae - then syls that start with lower FR would
                % show artificial increase in FR (and vice versa for those
                % starting higher) - quick solution, scale all
                % drift-subtracted values by mean FR
                
                ThisNSpksMat_minusDiff = ThisNSpksMat_minusDiff./repmat(ThisNspksDrift, 1, size(ThisNSpksMat,2));
                end
                
                corr_FF_vs_SpkMinusDiff = nan(1,nmotifs);
                for j=1:nmotifs
                    
                    ffvals = ThisFFmat(ThistrainIndsSongs,j);
                    nspk = ThisNSpksMat_minusDiff(ThistrainIndsSongs,j);
                    
                    rho = corr(nspk, ffvals);
                    
                    corr_FF_vs_SpkMinusDiff(j) = rho;
                end
                
                
                % ---------------- BASE CORR (FF vs. Nspks)
                basecorr_NspkminusdriftVsFF = nan(1,nmotifs);
                basecorr_NspkminusdriftVsFF_p = nan(1,nmotifs);
                for j=1:nmotifs
                    ffvals = ThisFFmat(ThisbaseIndsSongs,j);
                    nspks = ThisNSpksMat_minusDiff(ThisbaseIndsSongs, j);
                    
                    if all(isnan(ffvals))
                        continue
                    end
                    [rho, p] = corr(ffvals, nspks);
                    basecorr_NspkminusdriftVsFF(j) = rho;
                    basecorr_NspkminusdriftVsFF_p(j) = p;
                end
                
                AllBaseFFvsSpkminusdiff_corr = [AllBaseFFvsSpkminusdiff_corr; basecorr_NspkminusdriftVsFF'];
                AllBaseFFvsSpkminusdiff_corr_p = [AllBaseFFvsSpkminusdiff_corr_p; basecorr_NspkminusdriftVsFF_p'];
                
                
                % ------------------ Nspks (train End minus start)
                spkdiff = nan(1,nmotifs);
                spkdiff_frac = nan(1,nmotifs);
                ThistrainIndsSongs_Late = ThistrainIndsSongs(ceil(end/2):end);
                ThistrainIndsSongs_Early = ThistrainIndsSongs(1:ceil(end/2)-1);
                
                for j=1:nmotifs
                    
                    nspkL = ThisNSpksMat(ThistrainIndsSongs_Late,j);
                    nspkE = ThisNSpksMat(ThistrainIndsSongs_Early, j);
%                     nspkE = ThisNSpksMat(ThisbaseIndsSongs, j);
                    
                    spkdiff(j) = median(nspkL) - median(nspkE);
                    spkdiff_frac(j) = (median(nspkL) - median(nspkE))/(median(nspkE));
                    
                end
                
                
                % ------------------ Nspks (train End minus Base) (minus
                % Diff)
                spkdiff_TrainMinusBase_MinusDiff = nan(1,nmotifs);
                for j=1:nmotifs
                    
                    nspkL = ThisNSpksMat_minusDiff(ThistrainIndsSongs_Late,j);
                    nspkE = ThisNSpksMat_minusDiff(ThisbaseIndsSongs, j);
                    
                    spkdiff_TrainMinusBase_MinusDiff(j) = median(nspkL) - median(nspkE);
                end
                
                
                
                % ---------------------- xcorr of spikes (minus diff syl)
                % vs. FF
                xcorr_NspkminusDiff_vs_FF = {};
                for j=1:nmotifs
                    
                    ffvals = ThisFFmat(ThistrainIndsSongs,j);
                    nspk = ThisNSpksMat_minusDiff(ThistrainIndsSongs,j);
                    
                    [cc, lags] = xcov(nspk, ffvals, 10, 'Coeff');
                                        
                    xcorr_NspkminusDiff_vs_FF = [xcorr_NspkminusDiff_vs_FF cc'];
                end

                
                
                
               %% ================ OUTPUT
                AllIsTarg = [AllIsTarg; ThisIsTargVec'];
                AllIsSame = [AllIsSame; ThisIsSameVec'];
                AllBirdnum = [AllBirdnum; i*ones(size(ThisIsTargVec'))];
                AllExptnum = [AllExptnum; ii*ones(size(ThisIsTargVec'))];
                AllSwnum = [AllSwnum;  iii*ones(size(ThisIsTargVec'))];
                AllMotifnum = [AllMotifnum; ThisMotifOrigNums'];
                AllNeurnum = [AllNeurnum; nn*ones(size(ThisIsTargVec'))];
                AllTrainContingency = [AllTrainContingency; repmat(learnconting, length(ThisIsTargVec),1)];
                
                % ------------- RAW DAT
                tmp = mat2cell(ThisNSpksMat', ones(1, size(ThisNSpksMat,2)), size(ThisNSpksMat,1));
                AllRawNspks = [AllRawNspks; tmp];
                
                tmp = mat2cell(ThisNSpksMat_minusDiff', ones(1, size(ThisNSpksMat_minusDiff,2)), size(ThisNSpksMat_minusDiff,1));
                AllRawNspksMinusDrift = [AllRawNspksMinusDrift; tmp];

                tmp = ThisFRsmooth';
                for j=1:size(tmp,1)
                   AllRawFRSmooth = [AllRawFRSmooth; cell2mat(tmp(j,:))];
                end
                
                tmp = mat2cell(ThisFFmat', ones(1, size(ThisFFmat,2)), size(ThisFFmat,1));
                AllRawFF = [AllRawFF; tmp];
                
                tmp = mat2cell(repmat(ThisTvalsVec, size(ThisNSpksMat,2),1), ...
                    ones(1, size(ThisNSpksMat,2)), size(ThisTvalsVec,2));
                AllRawTvals = [AllRawTvals; tmp];
                
                tmp = mat2cell(repmat(ThistrainIndsSongs, size(ThisNSpksMat,2),1), ...
                    ones(1, size(ThisNSpksMat,2)), size(ThistrainIndsSongs,2));
                AllRawTrainInds = [AllRawTrainInds; tmp];
                
                 tmp = mat2cell(repmat(ThisbaseIndsSongs, size(ThisNSpksMat,2),1), ...
                    ones(1, size(ThisNSpksMat,2)), size(ThisbaseIndsSongs,2));
                AllRawBaseInds = [AllRawBaseInds; tmp];
                
                
               % ----
               
                AllFFslope = [AllFFslope; learnslope'];
                AllFFslopeSig = [AllFFslopeSig; learnsig'];
                AllSpkSTDsongs = [AllSpkSTDsongs; spkstd'];
                AllFFvsSpk_corr = [AllFFvsSpk_corr; corr_FFvsSpk'];
                AllSpkDiff_lateminusearly = [AllSpkDiff_lateminusearly; spkdiff'];
                AllSpkDiffFrac_lateminusearly = [AllSpkDiffFrac_lateminusearly; spkdiff_frac'];
                
                AllFFvsSpkminusdiff_corr = [AllFFvsSpkminusdiff_corr; corr_FF_vs_SpkMinusDiff'];
                
                AllSpkDiff_minusDiff_lateminusbase = [AllSpkDiff_minusDiff_lateminusbase; ...
                    spkdiff_TrainMinusBase_MinusDiff'];
                
                AllXcorr_NspkminusDiff_vs_FF = [AllXcorr_NspkminusDiff_vs_FF; xcorr_NspkminusDiff_vs_FF'];
                
              %% ================= get all pairwise correlations between targ and same
                
%                 corrAllPairs_FF_vs_SpkMinusDiff = nan(nmotifs,nmotifs); % ff of dim1 vs nspk of dim 2
                corrAllPairs_FF_vs_SpkMinusDiff = cell(nmotifs, 1);
                for j=1:nmotifs
                    corrAllPairs_FF_vs_SpkMinusDiff{j} = nan(1,nmotifs);
                    
                    for jj=1:nmotifs
                        
                    ffvals = ThisFFmat(ThistrainIndsSongs,j);
                    nspk = ThisNSpksMat_minusDiff(ThistrainIndsSongs,jj);
                    
                    rho = corr(nspk, ffvals);
                    
%                     corrAllPairs_FF_vs_SpkMinusDiff(j, jj) = rho;    
                    corrAllPairs_FF_vs_SpkMinusDiff{j}(jj) = rho;
                    end
                end
                
                AllPairsMotifs_FF_vs_SpkMinusDiff_corr = [AllPairsMotifs_FF_vs_SpkMinusDiff_corr; ...
                    corrAllPairs_FF_vs_SpkMinusDiff];
                
              
                
                %% ======= PLOT FOR THIS NEURON
                
                if strcmp(birdname, birdname_get) & ...
                        strcmp(exptname, exptname_get) & ...
                        switchnum_get == iii
                    
                    
                    figcount=1;
                    subplotrows=6;
                    subplotcols=1;
                    fignums_alreadyused=[];
                    hfigs=[];
                    hsplots = [];
                    
                    for j=1:size(ThisFFmat,2)
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(nn) '-mot' num2str(j)]);
                        ylabel('ff');
                        ff = ThisFFmat(:,j);
                        
                        % -- plotcol
                        if ThisIsTargVec(j)==1
                            pcol = 'k';
                        elseif ThisIsSameVec(j)==1
                            pcol = 'b';
                        else
                            pcol = 'r';
                        end
                        
                        if all(isnan(ff))
                            lt_plot_annotation(1, ['skipped motif' num2str(j) '(nan)'], 'r')
                        else
                            plot(ThisTvalsVec, ff, 'o', 'Color', pcol);
                        end
                        line([ThisTvalsVec(ThisbaseIndsSongs(end)) ThisTvalsVec(ThisbaseIndsSongs(end)) ], ...
                            ylim);
                    end
                    
                    axis tight
                    linkaxes(hsplots, 'xy');
                    
                    % -----------------
                    figcount=1;
                    subplotrows=6;
                    subplotcols=1;
                    fignums_alreadyused=[];
                    hfigs=[];
                    hsplots = [];
                    
                    for j=1:size(ThisNSpksMat,2)
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        hsplots = [hsplots hsplot];
                        title('nspks');
                        y = ThisNSpksMat(:,j);
                        
                        % -- plotcol
                        if ThisIsTargVec(j)==1
                            pcol = 'k';
                        elseif ThisIsSameVec(j)==1
                            pcol = 'b';
                        else
                            pcol = 'r';
                        end
                        
                        plot(ThisTvalsVec, y, 'o', 'Color', pcol);
                        line([ThisTvalsVec(ThisbaseIndsSongs(end)) ThisTvalsVec(ThisbaseIndsSongs(end)) ], ...
                            ylim);
                    end
                    axis tight
                    linkaxes(hsplots, 'xy');
                    
                end
                
                
                
            end
        end
    end
end

AllTrainDir = sign(AllTrainContingency(:, 2) - AllTrainContingency(:, 1));
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

%% ======================== PLOT ALL EXPERIMENTS [timecourses]

maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

plotTime =0; % if 0, tjhen rends
onlyPlotGoodLearn =0;
plotEachSylSep = 1; % if 1, then plots Nspk(minus drift) vs. FF
PlotSmFr=1; % if 1, then plots sm fr at end compare to baseline.

for i=5
%     for ii=1:maxexpts
        for ii=3
%         for iii=1:maxswitch
            for iii=1:maxswitch
            plotted = 0;
            for j=1:maxneuron
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                
                % ---- skip if not good learning
                if onlyPlotGoodLearn==1
                    inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1;
                    if all(AllFFslopeSig(inds)==0)
                        disp('SKIP - NOT GOOD LEARNING');
                        % --- then skip, since not good learning
                        continue
                    end
                end
                
                bname = MOTIFSTATS_Compiled.birds(i).birdname;
                exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
                hsplots = [];
                lt_figure; hold on;
                
                % ========================== OVERLAY LEARNING
                hsplot = lt_subplot(4,1,1); hold on;
                hsplots = [hsplots hsplot];
                ylabel('ff');
                title([bname '-' exptname '-sw' num2str(iii) '-n' num2str(j)]);
                
                % ---- TARG
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                plotcol = 'k';
                
                ffmat = cell2mat(AllRawFF(inds));
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspks(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                baseInds = cell2mat(AllRawBaseInds(inds));
                trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
                learn = AllFFslope(inds);
                learnsig = AllFFslopeSig(inds);
                
                for k=1:size(ffmat,1)
                    if plotTime==0
                        plot(ffmat(k,[baseInds trainInds]), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), ffmat(k,[baseInds trainInds]), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                    end
                    
                    % --- overlay if sig learning
                    if learnsig(k)==1
                        pcol = 'm';
                    else
                        pcol = 'k';
                    end
                    lt_plot_text(0, 15*k, ['ffslope:' num2str(learn(k))], pcol);
                end
                if plotTime==1
                    line([tvals(1, trainInds(1)-1) tvals(1, trainInds(1)-1)], ylim);
                else
%                     line([trainInds(1)-1 trainInds(1)-1], ylim);
                    line([length(baseInds) length(baseInds)], ylim);
                end
                
                % ---- SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                plotcol = 'b';
                
                ffmat = cell2mat(AllRawFF(inds));
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspks(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                baseInds = cell2mat(AllRawBaseInds(inds));
                trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                

                for k=1:size(ffmat,1)
                    if plotTime==0
                        plot(ffmat(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), ffmat(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    end
                end
                
                
                
                % ---- DIFF
inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==0;
                plotcol = 'r';
                
                ffmat = cell2mat(AllRawFF(inds));
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspks(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                baseInds = cell2mat(AllRawBaseInds(inds));
                trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
                for k=1:size(ffmat,1)
                    if plotTime==0
                        plot(ffmat(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), ffmat(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    end
                end
                
                lt_plot_zeroline;
                
                
                % =============================== NSPKS
                hsplot = lt_subplot(4,1,2); hold on;
                hsplots = [hsplots hsplot];
                ylabel('nspks');
                
                % ---- TARG
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                plotcol = 'k';
                
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspks(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                baseInds = cell2mat(AllRawBaseInds(inds));
                   trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
             
                for k=1:size(nspks,1)
                    if plotTime==0
                        plot(nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                    end
                end
                if plotTime==1
                    line([tvals(1, trainInds(1)-1) tvals(1, trainInds(1)-1)], ylim);
                else
                    line([length(baseInds) length(baseInds)], ylim);
%                     line([trainInds(1)-1 trainInds(1)-1], ylim);
                end
                
                % ---- SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                plotcol = 'b';
                
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspks(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                 baseInds = cell2mat(AllRawBaseInds(inds));
                   trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
            
                for k=1:size(nspks,1)
                    if plotTime==0
                        plot(nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    end
                end
                
                
                % ---- DIFF
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==0;
                plotcol = 'r';
                
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspks(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                 baseInds = cell2mat(AllRawBaseInds(inds));
                      trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
          
                for k=1:size(nspks,1)
                    if plotTime==0
                        plot(nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    end
                end
                
                lt_plot_zeroline;
                
                
                % ========================================= NSPKS, minus
                % drift
                hsplot = lt_subplot(4,1,3); hold on;
                hsplots = [hsplots hsplot];
                ylabel('nspks (minus drift)');
                
                % ---- TARG
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                plotcol = 'k';
                
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspksMinusDrift(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                 baseInds = cell2mat(AllRawBaseInds(inds));
                  trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
              
                for k=1:size(nspks,1)
                    if plotTime==0
                        plot(nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                    end
                end
                if plotTime==1
                    line([tvals(1, trainInds(1)-1) tvals(1, trainInds(1)-1)], ylim);
                else
                    line([length(baseInds) length(baseInds)], ylim);
%                     line([trainInds(1)-1 trainInds(1)-1], ylim);
                end
                
                % ---- SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                plotcol = 'b';
                
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspksMinusDrift(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                 baseInds = cell2mat(AllRawBaseInds(inds));
                     trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
           
                for k=1:size(nspks,1)
                    if plotTime==0
                        plot(nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    end
                end
                
                
                % ---- DIFF
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==0;
                plotcol = 'r';
                
                tvals = cell2mat(AllRawTvals(inds));
                nspks = cell2mat(AllRawNspksMinusDrift(inds));
                trainInds = cell2mat(AllRawTrainInds(inds));
                 baseInds = cell2mat(AllRawBaseInds(inds));
                  trainInds = trainInds(1,:);
                baseInds= baseInds(1,:);                
                
              
                for k=1:size(nspks,1)
                    if plotTime==0
                        plot(nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    else
                        plot(tvals(k,[baseInds trainInds]), nspks(k,[baseInds trainInds]), '-o', 'Color', plotcol);
                    end
                end
                
                lt_plot_zeroline;
                
                
                % ======================================
                linkaxes(hsplots, 'x');
                
                % ===================================== SUMMARY METRICS
                lt_subplot(4,2,7); hold on;
                title('nspk (minus drift) vs. ff');
                ylabel('xcorr');
                
                % -- targ
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                
                cc = cell2mat(AllXcorr_NspkminusDiff_vs_FF(inds));
                x = -floor(size(cc,2)/2):1:floor(size(cc,2)/2);
                plot(x, cc, '-k');
                
                % -- same
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                
                if any(inds)
                    cc = cell2mat(AllXcorr_NspkminusDiff_vs_FF(inds));
                    x = -floor(size(cc,2)/2):1:floor(size(cc,2)/2);
                    plot(x, cc, '-b');
                    
                    axis tight;
                    lt_plot_zeroline;
                    lt_plot_zeroline_vert;
                end
                
                
                % ==================================================
                if plotEachSylSep==1
                    figcount=1;
                    subplotrows=5;
                    subplotcols=1;
                    fignums_alreadyused=[];
                    hfigs=[];
                    
                    % -------------- TARG
                    inds = find(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    pcol = 'k';
                    
                    for jj=inds'
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title('targ');
                        ylabel('ff(fill); nspk minus drift (x2) (open)');
                        
                        ff = AllRawFF{jj};
                        nspk = AllRawNspksMinusDrift{jj};
                        trainInds = cell2mat(AllRawTrainInds(jj));
                        baseInds = cell2mat(AllRawBaseInds(jj));
                        
                        plot(ff([baseInds trainInds]), '-o', 'Color', pcol, 'MarkerFaceColor', pcol);
                        plot(2*nspk([baseInds trainInds]), '-o', 'Color', pcol);
                        line([trainInds(1)-1 trainInds(1)-1], ylim);
                        
                        % --- annotate correlations
                        lt_plot_annotation(1, ['corr = ' num2str(AllFFvsSpkminusdiff_corr(jj))], 'r');
                    end
                    
                    
                    % -------------- SAME
                    inds = find(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1);
                    pcol = 'b';
                    
                    for jj=inds'
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title('same');
                        ylabel('ff(fill); nspk minus drift (x2) (open)');
                        
                        ff = AllRawFF{jj};
                        nspk = AllRawNspksMinusDrift{jj};
                        trainInds = cell2mat(AllRawTrainInds(jj));
                        baseInds = cell2mat(AllRawBaseInds(jj));

                        plot(ff([baseInds trainInds]), '-o', 'Color', pcol, 'MarkerFaceColor', pcol);
                        plot(2*nspk([baseInds trainInds]), '-o', 'Color', pcol);
                        line([trainInds(1)-1 trainInds(1)-1], ylim);
                        
                        % --- annotate correlations
                        lt_plot_annotation(1, ['corr = ' num2str(AllFFvsSpkminusdiff_corr(jj))], 'r');
                        
                    end
                end
                
                
                % ========================== PLOT SMOOTHED FR
                if PlotSmFr==1
                    figcount=1;
                    subplotrows=4;
                    subplotcols=2;
                    fignums_alreadyused=[];
                    hfigs=[];
                    
                    % -------------- TARG
                    inds = find(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    pcol = 'k';
                    
                    for jj=inds'
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title('targ');
                        
                        % --- get mean smooth FR during and base and train
                        frsmooth = AllRawFRSmooth{jj};
                        baseinds = AllRawBaseInds{jj};
                        traininds = AllRawTrainInds{jj};
                        trainindsLate = traininds(ceil(end/2):end);
                        
                        frbase = mean(frsmooth(:, baseinds), 2);
                        frbase_sem = lt_sem(frsmooth(:, baseinds)');
                        frtrain = mean(frsmooth(:, trainindsLate), 2);
                        frtrain_sem = lt_sem(frsmooth(:, trainindsLate)');
                        x = premotorWind(1)+0.001:0.001:premotorWind(end);
                        
                        shadedErrorBar(x, frbase, frbase_sem, {'Color', 'k'}, 1);
                        shadedErrorBar(x, frtrain, frtrain_sem, {'Color', 'r'}, 1);
                        
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                    end
                    
                    % -------------- SAME
                    inds = find(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1);
                    pcol = 'k';
                    
                    for jj=inds'
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title('same');
                        
                        % --- get mean smooth FR during and base and train
                        frsmooth = AllRawFRSmooth{jj};
                        baseinds = AllRawBaseInds{jj};
                        traininds = AllRawTrainInds{jj};
                        trainindsLate = traininds(ceil(end/2):end);
                        
                        frbase = mean(frsmooth(:, baseinds), 2);
                        frbase_sem = lt_sem(frsmooth(:, baseinds)');
                        frtrain = mean(frsmooth(:, trainindsLate), 2);
                        frtrain_sem = lt_sem(frsmooth(:, trainindsLate)');
                        x = premotorWind(1)+0.001:0.001:premotorWind(end);
                        
                        shadedErrorBar(x, frbase, frbase_sem, {'Color', 'k'}, 1);
                        shadedErrorBar(x, frtrain, frtrain_sem, {'Color', 'r'}, 1);
                        
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                    end
                    
                    
                    
                end
                
                plotted=1;
            end
            
            if plotted==1
                if strcmp(input('close? (y, or n)', 's'), 'y')
                close all;
                end
            end
        end
    end
end


%% ========== [HETEROGENEITY] how heterogenous are metrics?

% onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =0;
plottype = 1;
% % onlyIfStartWNoff =0;
% overlayIfSU=1;
% CombineWithinSwitch = 1;
% takeAbs = 1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);
% 
% X= [];
% Y = [];
% % Birdnum = [];
% SwitchCount = [];
% swcount = 1;

% ================= PICK METRIC TO PLOT
Ymetric = [];
ylabstr = '';
YLIM = [];
switch plottype
    case 0
        % ff vs. spk(minus drift) correaltion
        ylabstr = 'corr (nspkminusdrift vs. ff)';
        Ymetric = AllFFvsSpkminusdiff_corr;
        YLIM = [-1 1];
    case 1
        % -- change in nspk (minus drift)
        ylabstr = 'change in nspk (minus drift)';
        Ymetric = AllSpkDiff_minusDiff_lateminusbase;
end


for i=5
    lt_figure; hold on ;
    subplot(4,1,1:2); hold on;
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    title(birdname);
    xlabel('neuron # (fake) (red lines switch; magenta line expt) (fill: learn sig)');
    ylabel(ylabstr);
    ncount = 1;
    for ii=1:maxexpts
        for iii=1:maxswitch
            %             swcount = swcount+1;
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                %                 if onlyIfStartWNoff==1
                %                     if unique(All_StartFromWNOff(inds))==0
                %                         continue
                %                     end
                %                 end
                %
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                %
                %                 % skip if neuron does not contribute both targ and nontarg
                %                 if length(unique(AllIsTarg(inds)))==1
                %                     continue
                %                 end
                %
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    %                     if ~any(learnsig==1 & learntarg>0)
                    if any(learnsig==0) | any(learntarg<0)
                        continue
                    end
                    %                     disp(learntarg);
                end
                
                % ================== DISPLAY USEFUL THINGS
                % -------------- direction of learning
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii;
                learnconting = unique(AllTrainContingency(inds, :), 'rows');
                assert(size(learnconting,1)==1, 'asdf');
                
                
                exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
                disp([birdname '-' exptname '-sw' num2str(iii) '-n' num2str(j) 'conting:  ' num2str(learnconting)]);
                
                %                 if strcmp(birdname, 'wh6pk36') & strcmp(exptname, 'LMANlearn2') & iii==1
                %                     keyboard
                %                 end
                
                % =============== TARGS
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                pcol = 'k';
                plot(ncount, Ymetric(inds), 'o', 'Color', pcol);
                % --- write motif num
                mnums = AllMotifnum(inds);
                learnsig = AllFFslopeSig(inds);
                learnlearn = AllFFslope(inds);
                
                 y = Ymetric(inds);
                for jj=1:length(mnums)
                    lt_plot_text(ncount, y(jj), ['m' num2str(mnums(jj))], pcol);
                end
                % --- note down learning
                subplot(4,1,3); hold on;
                ylabel('laerning');
                for jj=1:length(learnlearn)
                    if learnsig(jj)==1
                   plot(ncount, learnlearn(jj), 'o', 'Color', pcol, 'MarkerFaceColor', pcol);
                    else
                        plot(ncount, learnlearn(jj), 'o', 'Color', pcol);
                    end
                    lt_plot_text(ncount, learnlearn(jj), [num2str(learnconting)], 'r', 5);
                end
                % -- note down base corr
                basecorr = AllBaseFFvsSpkminusdiff_corr(inds);   
                basecorr_p = AllBaseFFvsSpkminusdiff_corr_p(inds);
                subplot(4,1,4); hold on;
                ylabel('base FF/nspk corr');
                for jj=1:length(basecorr)
                    if basecorr_p(jj)<0.05
                        plot(ncount, basecorr(jj), 'o', 'Color', pcol, 'MarkerFaceColor', pcol);
                    else
                        plot(ncount, basecorr(jj), 'o', 'Color', pcol);
                    end
                end
                
                
                subplot(4,1,1:2);
                
                % =============== SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                if any(inds)
                    pcol = 'b';
                    plot(ncount, Ymetric(inds), 'o', 'Color', pcol);
                    % --- write motif num
                learnsig = AllFFslopeSig(inds);
                learnlearn = AllFFslope(inds);
                mnums = AllMotifnum(inds);
                    
                    y = Ymetric(inds);
                    for jj=1:length(mnums)
                        lt_plot_text(ncount, y(jj), ['m' num2str(mnums(jj))], pcol);
                        if learnsig(jj)==1
                            plot(ncount, y(jj), 'o', 'MarkerFaceColor', 'm');
                        end
                    end
                
                                    % --- note down learning
                subplot(4,1,3); hold on;
                for jj=1:length(learnlearn)
                    if learnsig(jj)==1
                   plot(ncount, learnlearn(jj), 'o', 'Color', pcol, 'MarkerFaceColor', pcol);
                    else
                        plot(ncount, learnlearn(jj), 'o', 'Color', pcol);
                    end
                end
                % -- note down base corr (spk, ff)
                basecorr = AllBaseFFvsSpkminusdiff_corr(inds);   
                basecorr_p = AllBaseFFvsSpkminusdiff_corr_p(inds);
                subplot(4,1,4); hold on;
                ylabel('base FF/nspk corr');
                for jj=1:length(basecorr)
                    if basecorr_p(jj)<0.05
                        plot(ncount, basecorr(jj), 'o', 'Color', pcol, 'MarkerFaceColor', pcol);
                    else
                        plot(ncount, basecorr(jj), 'o', 'Color', pcol);
                    end
                end

                
                subplot(4,1,1:2);

                end
                
                ncount = ncount+1;
            end
            
            
            
            line([ncount-0.5 ncount-0.5], ylim, 'Color', 'r');
        end
        line([ncount-0.5 ncount-0.5], ylim, 'Color', 'm');
    end
    xlim([0 ncount+1]);
    if isempty(YLIM)
%        axis tight 
    else
    ylim(YLIM);
    end
    lt_plot_zeroline;
    subplot(4,1,3); hold on;
    xlim([0 ncount+1]);
    lt_plot_zeroline;
    subplot(4,1,4); hold on;
    xlim([0 ncount+1]);
    lt_plot_zeroline;
    
    %     pause
    %     close all;
end


%% ========== [MORE RAW] [DISTRIBUTIONS, FILTERED] filtered by 

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

AllTrainDir = sign(AllTrainContingency(:, 2) - AllTrainContingency(:, 1));

for i=1:maxbirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
        % ------------ collect for global plot
        X = [];
        Y = [];
        TrainDir = [];
        StartFromWNoff = [];
        EndOnWNoff = [];
        SylType =[]; % 0 for targ, 1 for same
        NuerNum = [];

    for ii=1:maxexpts
        
        for iii=1:maxswitch
            
            inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii;
            if ~any(inds)
                continue
            end
            
            exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
            % =========== 1) Train DIR up
            if any(AllBirdnum==i & AllExptnum==ii & AllSwnum==iii & ...
                    AllTrainDir==1)==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                xlabel('base FF-spk corr');
                ylabel('Nspks change (red=WNoff)');
                title([birdname '-' exptname '-sw' num2str(iii) '-Train UP']);
                
                for j=1:maxneuron
                    inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==1;
                    
                    x = AllBaseFFvsSpkminusdiff_corr(inds);
                    y = AllSpkDiff_minusDiff_lateminusbase(inds);
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    wnoffStart = AllTrainContingency(inds,1)==0;
                    if unique(wnoffStart)==1
                    plot(x,y, '-', 'Color', [0.8 0.2 0.2]);    
                    else
                    plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
                    end
                    
                    for k=1:length(x)
                        if istarg(k) ==1
                            plot(x(k), y(k),'o', 'Color', 'k');
                        elseif istarg(k)==0 & issame(k)==1
                            plot(x(k), y(k), 'o', 'Color', 'b');
                        else
                            assert(1==2, 'asdfasd');
                        end
                    end
                    
                    lt_plot_zeroline;
                    lt_plot_zeroline_vert;
                    
                    xlim([min(AllBaseFFvsSpkminusdiff_corr) max(AllBaseFFvsSpkminusdiff_corr)]);
                    ylim([min(AllSpkDiff_minusDiff_lateminusbase) max(AllSpkDiff_minusDiff_lateminusbase)]);
                    
                    % ---------------------- COLLECT
                    X = [X; x];
                    Y = [Y; y];
                    TrainDir = [TrainDir; 1*ones(length(x),1)];
                    StartFromWNoff = [StartFromWNoff; wnoffStart];
                    EndOnWNoff = [EndOnWNoff; AllTrainContingency(inds,2)==0];
                    SylType =[SylType; issame]; % 0 for targ, 1 for same
                    NuerNum = [NuerNum; AllNeurnum(inds)];
                end
            end
            
            % =========== 1) Train DIR DOWN
            if any(AllBirdnum==i & AllExptnum==ii & AllSwnum==iii & ...
                    AllTrainDir==1)==-1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                xlabel('base FF-spk corr');
                ylabel('Nspks change (red=WNoff)');
                title([birdname '-' exptname '-sw' num2str(iii) '-Train UP']);
                
                for j=1:maxneuron
                    inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==1;
                    
                    x = AllBaseFFvsSpkminusdiff_corr(inds);
                    y = AllSpkDiff_minusDiff_lateminusbase(inds);
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    wnoffStart = AllTrainContingency(inds,1)==0;
                    if unique(wnoffStart)==1
                    plot(x,y, '-', 'Color', [0.8 0.2 0.2]);    
                    else
                    plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
                    end
                    
                    for k=1:length(x)
                        if istarg(k) ==1
                            plot(x(k), y(k),'o', 'Color', 'k');
                        elseif istarg(k)==0 & issame(k)==1
                            plot(x(k), y(k), 'o', 'Color', 'b');
                        else
                            assert(1==2, 'asdfasd');
                        end
                    end
                    
                    lt_plot_zeroline;
                    lt_plot_zeroline_vert;
                    
                    xlim([min(AllBaseFFvsSpkminusdiff_corr) max(AllBaseFFvsSpkminusdiff_corr)]);
                    ylim([min(AllSpkDiff_minusDiff_lateminusbase) max(AllSpkDiff_minusDiff_lateminusbase)]);
                    
                    % ---------------------- COLLECT
                    X = [X; x];
                    Y = [Y; y];
                    TrainDir = [TrainDir; 1*ones(length(x),1)];
                    StartFromWNoff = [StartFromWNoff; wnoffStart];
                    EndOnWNoff = [EndOnWNoff; AllTrainContingency(inds,2)==0];
                    SylType =[SylType; issame]; % 0 for targ, 1 for same
                    NuerNum = [NuerNum; AllNeurnum(inds)];
                end
            end

        end
    end
end

linkaxes(hsplots, 'xy');


%% ========== [IMPORTANT] [SAME AS ABOVE, PLOT BY [DISTRIBUTIONS, FILTERED] filtered by 

figcount=1;
subplotrows=6;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

AllTrainDir = sign(AllTrainContingency(:, 2) - AllTrainContingency(:, 1));

plotBaseCorr=0; % if 1, then makes x axis baseline pitch corr
averageWithinSwitch = 1; % if 0, then plots all motifs/neurons, 1, then one datapoint for each switch

for i=1:maxbirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    
    % =========== 1) Train DIR up [FROM WN OFF]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    xlabel('base FF-spk corr');
    ylabel('Nspks change (red=WNoff)');
    title([birdname ' OFF->UP']);
    
    for ii=1:maxexpts
        for iii=1:maxswitch
            for j=1:maxneuron
                inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==1 ...
                    & AllTrainContingency(:,1)==0;
                
                fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
                    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)
            end
        end
    end
    
    % =========== 1) Train DIR DN [FROM WN OFF]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    xlabel('base FF-spk corr');
    ylabel('Nspks change (red=WNoff)');
    title([birdname ' OFF->DN']);
    
    for ii=1:maxexpts
        for iii=1:maxswitch
            for j=1:maxneuron
                inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==-1 ...
                    & AllTrainContingency(:,1)==0;
                fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
                    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)
                
            end
        end
    end
    
    % =========== 1) Train DIR UP [SWITCH]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    xlabel('base FF-spk corr');
    ylabel('Nspks change (red=WNoff)');
    title([birdname ' DN->UP']);
    
    for ii=1:maxexpts
        for iii=1:maxswitch
            for j=1:maxneuron
                inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==1 ...
                    & AllTrainContingency(:,1)~=0 & AllTrainContingency(:,2)~=0;
                fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
                    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)
                
            end
        end
    end
    
    % =========== 1) Train DIR DN [SWITCH]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    xlabel('base FF-spk corr');
    ylabel('Nspks change (red=WNoff)');
    title([birdname ' UP->DN']);
    
    for ii=1:maxexpts
        for iii=1:maxswitch
            for j=1:maxneuron
                inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==-1 ...
                    & AllTrainContingency(:,1)~=0 & AllTrainContingency(:,2)~=0;
                fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
                    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)
            end
        end
    end
    
    % =========== 1) Train DIR UP [to WN OFF]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    xlabel('base FF-spk corr');
    ylabel('Nspks change (red=WNoff)');
    title([birdname ' DN->OFF']);
    
    for ii=1:maxexpts
        for iii=1:maxswitch
            for j=1:maxneuron
                inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==1 ...
                    & AllTrainContingency(:,2)==0;
                fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
                    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)
                
            end
        end
    end
    
    % =========== 1) Train DIR UP [to WN OFF]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    xlabel('base FF-spk corr');
    ylabel('Nspks change (red=WNoff)');
    title([birdname ' UP->OFF']);
    
    for ii=1:maxexpts
        for iii=1:maxswitch
            for j=1:maxneuron
                inds = AllBirdnum==i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & (AllIsTarg==1 | AllIsSame==1) & AllTrainDir==-1 ...
                    & AllTrainContingency(:,2)==0;
                fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
                    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)
            end
        end
    end
    
end

linkaxes(hsplots, 'xy');


%% ========== [MOTIF PAIRS] corr of FF vs. Nspk (minus drift)


onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 2;
% onlyIfStartWNoff =0;
overlayIfSU=1;
CombineWithinSwitch = 1;
takeAbs = 1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

X= [];
Y = [];
% Birdnum = [];
SwitchCount = [];
swcount = 1;
for i=1:maxbirds
    lt_figure; hold on ;
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    title(birdname);
    xlabel('neuron # (fake) (dashed=dat vs. shuff; red=switch; magenta=expt)');
    ylabel('corr (nspkminusdrift vs. ff)');
ncount = 1;
    for ii=1:maxexpts
        for iii=1:maxswitch
            swcount = swcount+1; 
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                %                 if onlyIfStartWNoff==1
                %                     if unique(All_StartFromWNOff(inds))==0
                %                         continue
                %                     end
                %                 end
                %
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                %
                %                 % skip if neuron does not contribute both targ and nontarg
                %                 if length(unique(AllIsTarg(inds)))==1
                %                     continue
                %                 end
                %
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    %                     if ~any(learnsig==1 & learntarg>0)
                    if any(learnsig==0) | any(learntarg<0)
                        continue
                    end
                    %                     disp(learntarg);
                end
                
                exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
                disp([birdname '-' exptname '-sw' num2str(iii)]);
                
%                 if strcmp(birdname, 'wh6pk36') & strcmp(exptname, 'LMANlearn2') & iii==1
%                     keyboard
%                 end
                
                % =============== TARGS
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                pcol = 'k';
                plot(ncount, AllFFvsSpkminusdiff_corr(inds), 'o', 'Color', pcol);
                % --- write motif num
                mnums = AllMotifnum(inds);
                y = AllFFvsSpkminusdiff_corr(inds);
                for jj=1:length(mnums)
                   lt_plot_text(ncount, y(jj), ['m' num2str(mnums(jj))], pcol);
                end
                
                % =============== SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                if any(inds)
                pcol = 'b';
                plot(ncount, AllFFvsSpkminusdiff_corr(inds), 'o', 'Color', pcol);
                % --- write motif num
                mnums = AllMotifnum(inds);
                y = AllFFvsSpkminusdiff_corr(inds);
                for jj=1:length(mnums)
                   lt_plot_text(ncount, y(jj), ['m' num2str(mnums(jj))], pcol);
                end
                end
                
                
                
                
                % ################################### SHUFFLED (targ - same
                % FF and premotor mixed)
                
                % ================= LOAD ALL MOTIFS FOR THIS NEURON
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                cmat = cell2mat(AllPairsMotifs_FF_vs_SpkMinusDiff_corr(inds));
                istarg = AllIsTarg(inds);
                issame = AllIsSame(inds);
                motifnums = AllMotifnum(inds);
                
                indsTarg = find(istarg==1);
                indsSame = find(istarg==0 & issame==1);
                
                % ================= TARG FF (same-type premotor)
                pcol = 'k';
                
                call = [];
                mnum_FF = []; 
                for k = indsTarg'
                    for kk = indsSame'
                   
                        call = [call cmat(k, kk)];
                        mnum_FF = [mnum_FF motifnums(k)];
                        
                    end
                end
                
                if ~isempty(call)
                plot(ncount+0.5, call, 'o', 'Color', pcol);
                  for jj=1:length(mnum_FF)
                   lt_plot_text(ncount+0.5, call(jj), ['m' num2str(mnum_FF(jj))], pcol);
                  end
                end
                
                % ================= SAME-TYPE FF (targ premotor)
                pcol = 'b';
                
                call = [];
                mnum_FF = []; 
                for k = indsSame'
                    for kk = indsTarg'
                   
                        call = [call cmat(k, kk)];
                        mnum_FF = [mnum_FF motifnums(k)];
                        
                    end
                end
                
                if ~isempty(call)
                
                plot(ncount+0.5, call, 'o', 'Color', pcol);
                  for jj=1:length(mnum_FF)
                   lt_plot_text(ncount+0.5, call(jj), ['m' num2str(mnum_FF(jj))], pcol);
                  end
                end
                  
                line([ncount+0.8 ncount+0.8], ylim, 'LineStyle', '--');
                
                ncount = ncount+1;

            end
            line([ncount-0.5 ncount-0.5], ylim, 'Color', 'r');
        end
        line([ncount-0.5 ncount-0.5], ylim, 'Color', 'm');
    end
    xlim([0 ncount+1]);
    ylim([-1 1]);
    lt_plot_zeroline;
%     pause
%     close all;
end


%% ========== [MOTIF PAIRS] SUMMARY PLOT
% ############ 1) Collects data into vector (see below)
% ############ 2) Plots separate line for each neuron

lt_figure; hold on;
lt_subplot(3,1,1); hold on;
xlabel('TARG(dat)--SAME(dat)--TARGmotif(shuff)--SAMEmotif(shuff)');
% ====================================================
AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind = nan(size(AllBirdnum)); % for each motif, correlates its FF with
% the spikes from all other motifs of different kind(ie. targ--nontarg;
% ignores differen type);
% then takes mean of the absolute value across those pairs


% ==================================================
% onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
% plottype = 2;
% onlyIfStartWNoff =0;
overlayIfSU=1;
CombineWithinSwitch = 1;
takeAbs = 1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

YAll = [];
SwAll = [];
% Birdnum = [];
SwitchCount = [];
swcount = 1;
for i=1:maxbirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
ncount = 1;
    for ii=1:maxexpts
        for iii=1:maxswitch
            swcount = swcount+1; 
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                %                 if onlyIfStartWNoff==1
                %                     if unique(All_StartFromWNOff(inds))==0
                %                         continue
                %                     end
                %                 end
                %
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                %
                %                 % skip if neuron does not contribute both targ and nontarg
                %                 if length(unique(AllIsTarg(inds)))==1
                %                     continue
                %                 end
                %
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    %                     if ~any(learnsig==1 & learntarg>0)
                    if any(learnsig==0) | any(learntarg<0)
                        continue
                    end
                    %                     disp(learntarg);
                end
                
                exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
                disp([birdname '-' exptname '-sw' num2str(iii)]);
                
%                 if strcmp(birdname, 'wh6pk36') & strcmp(exptname, 'LMANlearn2') & iii==1
%                     keyboard
%                 end
                
                
                % ################################### SHUFFLED (targ - same
                % FF and premotor mixed)
                 
                % ================= LOAD ALL MOTIFS FOR THIS NEURON
                inds = find(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j);
                
                cmat = cell2mat(AllPairsMotifs_FF_vs_SpkMinusDiff_corr(inds));
                istarg = AllIsTarg(inds);
                issame = AllIsSame(inds);
                motifnums = AllMotifnum(inds);
                
                indsTarg = find(istarg==1);
                indsSame = find(istarg==0 & issame==1);
                indsDiff = find(issame ==0);
                
                Cout = nan(size(cmat,1),1); % will still back in here
                
                if isempty(indsTarg) | isempty(indsSame)
                    continue
                end
               
                % ================= TARG FF (same-type premotor)
                call = [];
                mnum_FF = []; 
                for k = indsTarg'
                    cthis = [];
                    for kk = indsSame'
                   
                        cthis= [cthis cmat(k, kk)];
%                         mnum_FF = [mnum_FF motifnums(k)];
                        
                    end
                    call = [call; cthis];
                end
                
                % --- take average of absolute value
                call = mean(abs(call),2);
                Cout(indsTarg) = call;

                
                
                % ================ SAME
                call = [];
                mnum_FF = []; 
                for k = indsSame'
                    cthis = [];
                    for kk = indsTarg'
                   
                        cthis= [cthis cmat(k, kk)];
%                         mnum_FF = [mnum_FF motifnums(k)];
                        
                    end
                    call = [call; cthis];
                end
                
                % --- take average of absolute value
                call = mean(abs(call),2);
                Cout(indsSame) = call;

                
                % --------------- SANITY CHECKE - targ and sametype fileld
                % in
                all(~isempty(Cout([indsTarg; indsSame])));
                all(isempty(Cout([indsDiff])));
                
                
                % =========================== STICK BACK INTO OVERAL VECTOR
                AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind(inds) = Cout;
                
               
                
              %% ################## PLOT FOR THIS NEURON
                Y = nan(1,4); % targ(real)-same(real)-targ(shuff)-same(shuff)
                
                % --- targ real
              inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                Y(1) = mean(abs(AllFFvsSpkminusdiff_corr(inds)));
                
                % --- same real
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                Y(2) = mean(abs(AllFFvsSpkminusdiff_corr(inds)));
                
                % --- targ shuff
              inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                Y(3) = mean(AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind(inds));
                
                % --- same shuff
              inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                Y(4) = mean(AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind(inds));
                
                plot(1:4, Y, '-ok');
                
                
                YAll = [YAll; Y];
                SwAll = [SwAll; swcount];

                
                
            end
        end
%         line([ncount-0.5 ncount-0.5], ylim, 'Color', 'm');
    end
    ylim([-1 1]);
    lt_plot_zeroline;
%     pause
%     close all;
end

% ---------------------------------------------
xlim([0 5]);
ylim([0 1]);


% -------------------------------------------------------
YAll_bysw = grpstats(YAll, SwAll);
lt_subplot(3,1,2); hold on;
title('same as above, but avg for each switch');
x = 1:4;
plot(x, YAll_bysw, '-ok');
xlim([0 5]);
lt_plot_zeroline;




%% [PLOT CROSS-CORRELATIONS] =================

onlySigLearning = 1;
onlyIfHasSameSyl =1;
takeAbs=1;

maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);


CCall = cell(1,2);


for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            plotted = 0;
            CCswitch = [];
            for j=1:maxneuron
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                
                % skip if doesn't have same syl
                if onlyIfHasSameSyl==1
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                if ~any(inds)
                    continue
                end
                end
                
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    %                     if ~any(learnsig==1 & learntarg>0)
                    if ~any(learnsig==1 & learntarg>0)
                        continue
                    end
                    %                     disp(learntarg);
                end
                
                
                bname = MOTIFSTATS_Compiled.birds(i).birdname;
                exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
                hsplots = [];
                
                
                % ================================================
                % --------------- TARG
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                
                cc = cell2mat(AllXcorr_NspkminusDiff_vs_FF(inds));
                CCall{1} = [CCall{1}; cc];
                
                % --------------- SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                
                cc = cell2mat(AllXcorr_NspkminusDiff_vs_FF(inds));
                CCall{2} = [CCall{2}; cc];
                
                
            end
            
        end
    end
end

x = -floor(size(CCall{1},2)/2):1:floor(size(CCall{1},2)/2);
lt_figure; hold on; 

% --- targ
lt_subplot(3,1,1); hold on;
title('targ');
ind = 1;

plot(x, CCall{ind}, '-', 'Color', [0.7 0.7 0.7]);
plot(x, mean(CCall{ind},1), 'Color', 'k');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% --- targ
lt_subplot(3,1,2); hold on;
title('same');
ind = 2;

plot(x, CCall{ind}, '-', 'Color', [0.7 0.7 0.7]);
plot(x, mean(CCall{ind},1), 'Color', 'k');
lt_plot_zeroline;
lt_plot_zeroline_vert;

%% [COMPUTE NEW] ================= cross correlations between target nontargets
onlySigLearning = 1;

maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);


CCall = [];
LagsAll = [];

for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            plotted = 0;
            CCswitch = [];
            for j=1:maxneuron
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    %                     if ~any(learnsig==1 & learntarg>0)
                    if any(learnsig==0) | any(learntarg<0)
                        continue
                    end
                    %                     disp(learntarg);
                end
                
                
                bname = MOTIFSTATS_Compiled.birds(i).birdname;
                exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
                hsplots = [];
                %                 lt_figure; hold on;
                
                % ========================== OVERLAY LEARNING
                %                 hsplot = lt_subplot(2,1,1); hold on;
                %                 hsplots = [hsplots hsplot];
                %                 ylabel('ff');
                %                 title([bname '-' exptname '-sw' num2str(iii) '-n' num2str(j)]);
                
                % ---- TARG VS SAME TYPE
                inds1 = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                inds2 = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                
                ffmat1 = cell2mat(AllRawFF(inds1));
                ffmat2 = cell2mat(AllRawFF(inds2));
                nspk1 = cell2mat(AllRawNspks(inds1));
                nspk2 = cell2mat(AllRawNspks(inds2));
                %                 ffmat1 = nspk1;
                %                 ffmat2 = nspk2;
                trainInds = unique(cell2mat(AllRawTrainInds(inds)));
                baseInds = unique(cell2mat(AllRawBaseInds(inds)));
                
                % --- REMOVE DATAPOINTS THAT HAVE NAN'
                ffmat1(any(isnan(ffmat1)'), :) = [];
                ffmat2(any(isnan(ffmat2)'), :) = [];
                
                n1 = size(ffmat1,1);
                n2 = size(ffmat2, 1);
                %                 CCall = [];
                %                 LagsAll = [];
                
                for k=1:n1
                    for kk=1:n2
                        
                        % ============ training
                        [cc2, lags] = xcov(ffmat1(k,trainInds), ffmat2(kk,trainInds), 50, 'Coeff');
                        %                         CCall = [CCall; cc];
                        %                         LagsAll = [LagsAll; lags];
                        
                        % ============ baseline
                        [cc1, lags] = xcov(ffmat1(k,baseInds), ffmat2(kk,baseInds), 50, 'Coeff');
                        
                        cc = cc2-cc1;
                        CCswitch = [CCswitch; cc];
                    end
                end
                
                
                %                 for k=1:size(ffmat,1)
                %                     if plotTime==0
                %                    plot(ffmat(k,:), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                %                     else
                %                    plot(tvals(k,:), ffmat(k,:), '-o', 'Color', plotcol, 'MarkerFaceColor', plotcol);
                %                     end
                %
                %                     % --- overlay if sig learning
                %                     if learnsig(k)==1
                %                         pcol = 'm';
                %                     else
                %                         pcol = 'k';
                %                     end
                %                     lt_plot_text(0, 15*k, ['ffslope:' num2str(learn(k))], pcol);
                %                 end
                %                 if plotTime==1
                %                     line([tvals(1, trainInds(1)-1) tvals(1, trainInds(1)-1)], ylim);
                %                 else
                %                     line([trainInds(1)-1 trainInds(1)-1], ylim);
                %                 end
                
            end
            
            CCall = [CCall; mean(CCswitch,1)];
        end
    end
end

figure; hold on; plot(lags, CCall, '-', 'Color', [0.7 0.7 0.7]);
plot(lags, mean(CCall, 1), '-r')
lt_plot_zeroline;
lt_plot_zeroline_vert;

%% ========================= PLOT
% nontarg)
% onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 4;
% onlyIfStartWNoff =0;
% overlayIfSU=1;
CombineWithinSwitch = 0;
takeAbs = 0;
flipIfLearnDirNeg = 0;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

lt_figure; hold on;
lt_subplot(4, 1, 1:2); hold on;
title('TARG (k), sametype(b), diff(r) [each neuron paired]');
xlabel('learn rate');

X = [];
Y = [];
% Birdnum = [];
SwitchCount = [];
swcount = 1;
BirdCount = [];
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            swcount = swcount+1; 
            for j = 1:maxneuron
                
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j;
                
                %                 if onlyIfStartWNoff==1
                %                     if unique(All_StartFromWNOff(inds))==0
                %                         continue
                %                     end
                %                 end
                %
                
                % skip if no dat
                if ~any(inds)
                    continue
                end
                %
                %                 % skip if neuron does not contribute both targ and nontarg
                %                 if length(unique(AllIsTarg(inds)))==1
                %                     continue
                %                 end
                %
                % skip if no target has significant learning
                if onlySigLearning==1
                    learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                        & AllNeurnum==j & AllIsTarg==1);
                    
                    %                     if ~any(learnsig==1 & learntarg>0)
                    if any(learnsig==0) | any(learntarg<0)
                        continue
                    end
                    %                     disp(learntarg);
                end
                
                %                 if onlyKeepIfGreaterLearnNontarg==1
                %                     targlearn = mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                %                         & AllNeurnum==j & AllIsTarg==1));
                %                     nontarglearn =  mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                %                         & AllNeurnum==j & AllIsTarg==0));
                %
                %                     if targlearn>nontarglearn
                %                         continue
                %                     end
                %                 end
                
                
                % ----------- TO COLLECT
                thisLearn = [];
                thisYval = [];
                thisTarg = [];
                thisSame = [];
                
                % ======================================= TARGET
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==1;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = AllFFvsSpk_corr(inds);
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                        case 2
                            % 
                            yval = AllFFvsSpkminusdiff_corr(inds);
                            ylabel('FF vs SpkMinusDiffSyl, corr');
                        case 3
                            % 
                            yval = AllSpkDiff_minusDiff_lateminusbase(inds);
                            ylabel('late minus base (Nspk, minus diff syls)');
                        case 4
                            % 
                            yval = AllSpkDiff_lateminusearly(inds);
                            ylabel('late minus early (Nspk)');
                        case 5 
                            yval = AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind(inds);
                            ylabel('FF vs SpkMinusDiffSyl, corr (SHUFFLE spk vs. ff)');
                    end
                    % -- flip if learn dir neg
                    if flipIfLearnDirNeg==1
learndir = AllTrainDir(inds);
                    yval = yval.*learndir;
                end
                
                    if takeAbs==1
                        yval = abs(yval);
                    end
                    
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                    
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];
                
                
                % ======================================= SAME
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==1;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = AllFFvsSpk_corr(inds);
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                        case 2
                            % 
                            yval = AllFFvsSpkminusdiff_corr(inds);
                            ylabel('FF vs SpkMinusDiffSyl, corr');
                        case 3
                            % 
                            yval = AllSpkDiff_minusDiff_lateminusbase(inds);
                        case 4
                            % 
                            yval = AllSpkDiff_lateminusearly(inds);
                        case 5 
                            yval = AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind(inds);
                            ylabel('FF vs SpkMinusDiffSyl, corr (SHUFFLE spk vs. ff)');
                    end
                     % -- flip if learn dir neg
                    if flipIfLearnDirNeg==1
learndir = AllTrainDir(inds);
                    yval = yval.*learndir;
                end
                   
                    if takeAbs==1
                        yval = abs(yval);
                    end
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];
                
                
                
                
                % ======================================= DIFF
                inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllNeurnum==j & AllIsTarg==0 & AllIsSame==0;
                
                if ~any(inds)
                    % -- then give all nan
                    istarg = nan;
                    issame = nan;
                    learnrate = nan;
                    yval = nan;
                else
                    istarg = AllIsTarg(inds);
                    issame = AllIsSame(inds);
                    learnrate = AllFFslope(inds);
                    yval = nan;
                    switch plottype
                        case 0
                            % trial - spk vs, ff, corr
                            yval = AllFFvsSpk_corr(inds);
                            ylabel('spk vs. ff corr (abs value)');
                        case 1
                            % std of mean spks over trials
                            yval = AllSpkSTDsongs(inds);
                            ylabel('std of mean spks over songs (CV)');
                         case 2
                            % 
                            yval = AllFFvsSpkminusdiff_corr(inds);
                            ylabel('FF vs SpkMinusDiffSyl, corr');
                        case 3
                            % 
                            yval = AllSpkDiff_minusDiff_lateminusbase(inds);
                        case 4
                            % 
                            yval = AllSpkDiff_lateminusearly(inds);
                        case 5 
                            yval = AllPairsMotifs_FF_vs_SpkMinusDiff_corr_MeanAbsAcrossPreWind(inds);
                            ylabel('FF vs SpkMinusDiffSyl, corr (SHUFFLE spk vs. ff)');
                   end
                    % -- flip if learn dir neg
                    if flipIfLearnDirNeg==1
learndir = AllTrainDir(inds);
                    yval = yval.*learndir;
                end
                     if takeAbs==1
                        yval = abs(yval);
                    end
                   
                    learnrate = mean(learnrate(~isnan(learnrate)));
                    yval =  mean(yval(~isnan(learnrate)));
                    istarg = unique(istarg);
                    issame = unique(issame);
                end
                
                thisLearn = [thisLearn learnrate];
                thisYval = [thisYval yval];
                thisTarg = [thisTarg istarg];
                thisSame = [thisSame issame];
                
                
                
                
                
                X = [X; thisLearn];
                Y = [Y; thisYval];
%                 TargStat = [TargStat; thisTarg];
%                 SameStat = [SameStat; thisSame];
                SwitchCount = [SwitchCount; swcount];
                BirdCount = [BirdCount; i];
                
                % #################### count switch
                %                 SwitchCount = [SwitchCount swcounter];
                %                 swcounter = swcounter+1;
                %
                %                 Birdnum = [Birdnum i];
                
                % --------------- plot line connecting things
                %
                %                 if X(end)>X(end-1)
                %                     plotcol = 'b';
                %                 else
                %                     plotcol = [0.7 0.7 0.7];
                %                 end
                %                 if Y(end)<Y(end-1)
                %                     lstyle = '-';
                %                 else
                %                     lstyle = ':';
                %                 end
                %                 plot(X(end-1:end), Y(end-1:end), '-', 'Color', plotcol, ...
                %                     'LineStyle', lstyle);
                %
                
                %                 % ---- if plot SU
                %                 if overlayIfSU==1
                %                     if unique(All_SingleUnit(inds))==1
                %                         plot(X(end-1), Y(end-1), 'sm', 'MarkerSize', 15);
                %                         plot(X(end), Y(end), 'sm');
                %                     end
                %                 end
            end
        end
    end
end


% =============== if combine to get one datapt per switch
if CombineWithinSwitch==1
    Xnew = [];
    Ynew = [];
    Birdnew = [];
    enums = unique(SwitchCount);
    for j=enums'
        indtmp = SwitchCount==j;
        
        x = mean(X(indtmp,:),1);
        y = mean(Y(indtmp,:),1);
        bb = unique(BirdCount(indtmp));
        
        Xnew = [Xnew; x];
        Ynew = [Ynew; y];
        Birdnew = [Birdnew; bb];
    end
    X = Xnew;
    Y = Ynew;
    BirdCount = Birdnew;
end

% ================================ PLOT 1 (Y VS. LEARNING)
% lt_figure; hold on;
% lt_subplot(3, 1, 1:2); hold on;
title('TARG (r), nontarg(k) [each neuron paired]');
xlabel('learn rate');

% --------- FIRST: those that have paired data (targ vs. nontarg)
% indstmp = find(~any(isnan(X)'));
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% ----------- SECOND, plot those that are not paired (open circles)
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    %     plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    plot(X(j, 1), Y(j,1), 'ko');
    plot(X(j, 2), Y(j,2), 'bo');
    plot(X(j, 3), Y(j,3), 'ro');
end


YLIM = ylim;


% =============================================== PAIRED COMPARISON
lt_subplot(4,2,5); hold on;
xlabel('TARG -- SAMESYL -- DIFFSYL');
% ylabel('spk vs. ff corr (abs value)');
x = [1 2 3];
yy = [];

% --------- FIRST: those that have paired data (targ vs. nontarg)
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    plot(x, Y(j,:), 'o-k');
    %     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    %     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    %     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% --------- SECOND: not paired
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    plot(x+0.1, Y(j,:), 'o', 'Color', [0.7 0.7 0.7]);
    %     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    %     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    %     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

xlim([0 4]);

% --- sign rank
p = signrank(Y(:,1), Y(:,2));
lt_plot_text(1.5, 1.1*max(Y(:,1)), ['(1vs2) p=' num2str(p)]);
if ~all(isnan(Y(:,3)))
p = signrank(Y(:,1), Y(:,3));
lt_plot_text(2, 1.2*max(Y(:,1)), ['(1vs3) p=' num2str(p)]);
end



% =============================================== PAIRED COMPARISON
lt_subplot(4,2,6); hold on;
xlabel('learn diff (targ - nontarg)');
ylabel('diff in metric (targ - nontarg)');

y = Y(:,1) - Y(:,2);
x = X(:,1) - X(:,2);

lt_regress(y, x, 1, 0, 1, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ======================================= PAIRED, SEPARATED BY BIRD
lt_subplot(4,2, 7:8); hold on;
xlabel('TARG -- SAMESYL -- DIFFSYL (sep by bird)');

for i=1:maxbirds
    % ylabel('spk vs. ff corr (abs value)');
    x = 4*(i-1)+[1 2 3];
    
    % --------- FIRST: those that have paired data (targ vs. nontarg)
    indstmp = find(~any(isnan(X(:,1:2))') & BirdCount'==i); % pair those with both
    
    for j=indstmp
        plot(x, Y(j,:), 'o-k');
        %     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
        %     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
        %     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
    end
    
    % --------- SECOND: not paired
    indstmp = find(any(isnan(X(:,1:2))') & BirdCount'==i); % pair those with both
    
    for j=indstmp
        plot(x+0.1, Y(j,:), 'o', 'Color', [0.7 0.7 0.7]);
        %     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
        %     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
        %     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
    end
%     
%     xlim([0 4]);
%     
%     % --- sign rank
%     p = signrank(Y(:,1), Y(:,2));
%     lt_plot_text(1.5, 1.1*max(Y(:,1)), ['(1vs2) p=' num2str(p)]);
%     if ~all(isnan(Y(:,3)))
%         p = signrank(Y(:,1), Y(:,3));
%         lt_plot_text(2, 1.2*max(Y(:,1)), ['(1vs3) p=' num2str(p)]);
%     end
    
    lt_plot_zeroline;
end
%% ========================= PLOT [v2 -- SWITCH AS A DATPOINT]
% nontarg)
onlyKeepIfGreaterLearnNontarg=0;
onlySigLearning =1;
plottype = 1;
% onlyIfStartWNoff =0;
overlayIfSU=1;

% =================== first collect data so is matched - i.e. each neuron
% must contribute to both target and nontarget
maxbirds = max(AllBirdnum);
maxexpts = max(AllExptnum);
maxswitch = max(AllSwnum);
maxneuron = max(AllNeurnum);

lt_figure; hold on;
lt_subplot(3, 1, 1:2); hold on;
title('TARG (k), sametype(b), diff(r) [each neuron paired]');
xlabel('learn rate');

X = [];
Y = [];
% Birdnum = [];
TargStat = [];
SameStat = [];
% SwitchCount = [];
% swcounter = 1;
for i=1:maxbirds
    for ii=1:maxexpts
        for iii=1:maxswitch
            
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii;
            
            %                 if onlyIfStartWNoff==1
            %                     if unique(All_StartFromWNOff(inds))==0
            %                         continue
            %                     end
            %                 end
            %
            
            % skip if no dat
            if ~any(inds)
                continue
            end
            %
            %                 % skip if neuron does not contribute both targ and nontarg
            %                 if length(unique(AllIsTarg(inds)))==1
            %                     continue
            %                 end
            %
            % skip if no target has significant learning
            if onlySigLearning==1
                learnsig = AllFFslopeSig(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1);
                learntarg = AllFFslope(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                    & AllIsTarg==1);
                
                if ~any(learnsig==1 & learntarg>0)
                    continue
                end
                %                     disp(learntarg);
            end
            
            %                 if onlyKeepIfGreaterLearnNontarg==1
            %                     targlearn = mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
            %                         & AllNeurnum==j & AllIsTarg==1));
            %                     nontarglearn =  mean(AllLearnSlopeScaled(AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
            %                         & AllNeurnum==j & AllIsTarg==0));
            %
            %                     if targlearn>nontarglearn
            %                         continue
            %                     end
            %                 end
            
            
            % ----------- TO COLLECT
            thisLearn = [];
            thisYval = [];
            thisTarg = [];
            thisSame = [];
            
            % ======================================= TARGET
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii & AllIsTarg==1;
            
            if ~any(inds)
                % -- then give all nan
                istarg = nan;
                issame = nan;
                learnrate = nan;
                yval = nan;
            else
                istarg = AllIsTarg(inds);
                issame = AllIsSame(inds);
                learnrate = AllFFslope(inds);
                yval = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        yval = abs(AllFFvsSpk_corr(inds));
                        ylabel('spk vs. ff corr (abs value)');
                    case 1
                        % std of mean spks over trials
                        yval = AllSpkSTDsongs(inds);
                        ylabel('std of mean spks over songs (CV)');
                    case 2
                        yval = AllSpkDiff_lateminusearly(inds);
                        ylabel('diff in median spks (late minus early)');
                    case 3
                        yval = AllSpkDiffFrac_lateminusearly(inds);
                        ylabel('frac decrease in median spks (late minus early)');
                end
                
                learnrate = mean(learnrate(~isnan(learnrate)));
                yval =  mean(yval(~isnan(learnrate)));
                istarg = unique(istarg);
                issame = unique(issame);
            end
            
            thisLearn = [thisLearn learnrate];
            thisYval = [thisYval yval];
            thisTarg = [thisTarg istarg];
            thisSame = [thisSame issame];
            
            
            % ======================================= SAME
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==0 & AllIsSame==1;
            
            if ~any(inds)
                % -- then give all nan
                istarg = nan;
                issame = nan;
                learnrate = nan;
                yval = nan;
            else
                istarg = AllIsTarg(inds);
                issame = AllIsSame(inds);
                learnrate = AllFFslope(inds);
                yval = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        yval = abs(AllFFvsSpk_corr(inds));
                        ylabel('spk vs. ff corr (abs value)');
                    case 1
                        % std of mean spks over trials
                        yval = AllSpkSTDsongs(inds);
                        ylabel('std of mean spks over songs (CV)');
                    case 2
                        yval = AllSpkDiff_lateminusearly(inds);
                        ylabel('diff in median spks (late minus early)');
                    case 3
                        yval = AllSpkDiffFrac_lateminusearly(inds);
                        ylabel('frac decrease in median spks (late minus early)');
                end
                
                learnrate = mean(learnrate(~isnan(learnrate)));
                yval =  mean(yval(~isnan(learnrate)));
                istarg = unique(istarg);
                issame = unique(issame);
            end
            
            thisLearn = [thisLearn learnrate];
            thisYval = [thisYval yval];
            thisTarg = [thisTarg istarg];
            thisSame = [thisSame issame];
            
            
            
            
            % ======================================= DIFF
            inds = AllBirdnum == i & AllExptnum==ii & AllSwnum==iii ...
                & AllIsTarg==0 & AllIsSame==0;
            
            if ~any(inds)
                % -- then give all nan
                istarg = nan;
                issame = nan;
                learnrate = nan;
                yval = nan;
            else
                istarg = AllIsTarg(inds);
                issame = AllIsSame(inds);
                learnrate = AllFFslope(inds);
                yval = nan;
                switch plottype
                    case 0
                        % trial - spk vs, ff, corr
                        yval = abs(AllFFvsSpk_corr(inds));
                        ylabel('spk vs. ff corr (abs value)');
                    case 1
                        % std of mean spks over trials
                        yval = AllSpkSTDsongs(inds);
                        ylabel('std of mean spks over songs (CV)');
                    case 2
                        yval = AllSpkDiff_lateminusearly(inds);
                        ylabel('diff in median spks (late minus early)');
                    case 3
                        yval = AllSpkDiffFrac_lateminusearly(inds);
                        ylabel('frac decrease in median spks (late minus early)');
                end
                
                learnrate = mean(learnrate(~isnan(learnrate)));
                yval =  mean(yval(~isnan(learnrate)));
                istarg = unique(istarg);
                issame = unique(issame);
            end
            
            thisLearn = [thisLearn learnrate];
            thisYval = [thisYval yval];
            thisTarg = [thisTarg istarg];
            thisSame = [thisSame issame];
            
            
            
            
            
            X = [X; thisLearn];
            Y = [Y; thisYval];
            TargStat = [TargStat; thisTarg];
            SameStat = [SameStat; thisSame];
            
            
            
            % #################### count switch
            %                 SwitchCount = [SwitchCount swcounter];
            %                 swcounter = swcounter+1;
            %
            %                 Birdnum = [Birdnum i];
            
            % --------------- plot line connecting things
            %
            %                 if X(end)>X(end-1)
            %                     plotcol = 'b';
            %                 else
            %                     plotcol = [0.7 0.7 0.7];
            %                 end
            %                 if Y(end)<Y(end-1)
            %                     lstyle = '-';
            %                 else
            %                     lstyle = ':';
            %                 end
            %                 plot(X(end-1:end), Y(end-1:end), '-', 'Color', plotcol, ...
            %                     'LineStyle', lstyle);
            %
            
            %                 % ---- if plot SU
            %                 if overlayIfSU==1
            %                     if unique(All_SingleUnit(inds))==1
            %                         plot(X(end-1), Y(end-1), 'sm', 'MarkerSize', 15);
            %                         plot(X(end), Y(end), 'sm');
            %                     end
            %                 end
            
        end
    end
end

% ================================ PLOT 1 (Y VS. LEARNING)
% lt_figure; hold on;
% lt_subplot(3, 1, 1:2); hold on;
title('TARG (r), nontarg(k) [each neuron paired]');
xlabel('learn rate');

% --------- FIRST: those that have paired data (targ vs. nontarg)
% indstmp = find(~any(isnan(X)'));
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% ----------- SECOND, plot those that are not paired (open circles)
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    %     plot(X(j, :), Y(j,:), '-', 'Color', [0.7 0.7 0.7]);
    plot(X(j, 1), Y(j,1), 'ko');
    plot(X(j, 2), Y(j,2), 'bo');
    plot(X(j, 3), Y(j,3), 'ro');
end


YLIM = ylim;


% =============================================== PAIRED COMPARISON
lt_subplot(3,2,5); hold on;
xlabel('TARG -- SAMESYL -- DIFFSYL');
% ylabel('spk vs. ff corr (abs value)');
x = [1 2 3];
yy = [];

% --------- FIRST: those that have paired data (targ vs. nontarg)
indstmp = find(~any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    plot(x, Y(j,:), 'o-k');
    %     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    %     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    %     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

% --------- SECOND: not paired
indstmp = find(any(isnan(X(:,1:2))')); % pair those with both

for j=indstmp
    plot(x+0.1, Y(j,:), 'o', 'Color', [0.7 0.7 0.7]);
    %     lt_plot(X(j, 1), Y(j,1), {'Color', 'k'});
    %     lt_plot(X(j, 2), Y(j,2), {'Color', 'b'});
    %     lt_plot(X(j, 3), Y(j,3), {'Color', 'r'});
end

xlim([0 4]);

% --- sign rank
p = signrank(Y(:,1), Y(:,2));
lt_plot_text(1.5, 1.1*max(Y(:,1)), ['(1vs2) p=' num2str(p)]);
p = signrank(Y(:,1), Y(:,3));
lt_plot_text(2, 1.2*max(Y(:,1)), ['(1vs3) p=' num2str(p)]);



% =============================================== PAIRED COMPARISON
lt_subplot(3,2,6); hold on;
xlabel('learn diff (targ - nontarg)');
ylabel('diff in metric (targ - nontarg)');

y = Y(:,1) - Y(:,2);
x = X(:,1) - X(:,2);

lt_regress(y, x, 1, 0, 1, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;

%% FIT LME -
% NOTE: neurons numbers are specific to each experiment
% TO DO: combine multiple motifs for neuron - only keep neurons with both
% same and diff.

% ====== index each unique switch/expt/bird separately
dattable = table(AllBirdnum, AllExptnum, AllSwnum);
[C, ~, AllUniqueSwInd] = unique(dattable);


% ====== 
% inds = (AllIsSame==1 | AllIsTarg==1) & AllFFslopeSig==1; % only care about same and targ
inds = (AllIsSame==1 | AllIsTarg==1); % only care about same and targ
ds = table(AllSpkDiff_minusDiff_lateminusbase(inds), AllBirdnum(inds), AllExptnum(inds), ...
    AllSwnum(inds), AllNeurnum(inds), AllUniqueSwInd(inds), AllIsSame(inds), AllTrainDir(inds), ...
    'VariableNames',{'AllSpkDiff_minusDiff_lateminusbase', 'AllBirdnum', 'AllExptnum', ...
    'AllSwnum', 'AllNeurnum', 'AllUniqueSwInd', 'AllIsSame', 'AllTrainDir'});
formula = 'AllSpkDiff_minusDiff_lateminusbase ~ AllIsSame + (AllIsSame|AllUniqueSwInd)';
% formula = 'AllSpkDiff_minusDiff_lateminusbase ~ AllIsSame + (AllIsSame|AllUniqueSwInd) + (AllIsSame|AllUniqueSwInd:AllNeurnum)';
% formula = 'AllSpkDiff_minusDiff_lateminusbase ~ AllIsSame + AllIsSame:AllTrainDir + (AllIsSame|AllUniqueSwInd)';

lme = fitlme(ds, formula);





end


function fn_plotpair(AllBaseFFvsSpkminusdiff_corr, AllSpkDiff_minusDiff_lateminusbase, ...
    AllIsTarg, AllIsSame, inds, plotBaseCorr, averageWithinSwitch)


                if ~any(inds)
                    return
                end
                
                x = AllBaseFFvsSpkminusdiff_corr(inds);
                y = AllSpkDiff_minusDiff_lateminusbase(inds);
                istarg = AllIsTarg(inds);
                issame = AllIsSame(inds);
                
                if averageWithinSwitch==1
                   % -- separately take average
                   x = grpstats(x, issame);
                   y = grpstats(y, issame);
                   istarg = grpstats(istarg, issame);
                   issame = unique(issame);
                end
                
                
                if plotBaseCorr==1
                    plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
                    
                    for k=1:length(x)
                        if istarg(k) ==1
                            plot(x(k), y(k),'o', 'Color', 'k');
                        elseif istarg(k)==0 & issame(k)==1
                            plot(x(k), y(k), 'o', 'Color', 'b');
                        end
                    end
                    xlim([min(AllBaseFFvsSpkminusdiff_corr) max(AllBaseFFvsSpkminusdiff_corr)]);
                lt_plot_zeroline_vert;
                else
                    x = issame;
                    plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
                    plot(0, y(issame==0), 'ok');
                    if any(issame==1)
                    plot(1, y(issame==1), 'ob');
                    end
                    xlim([-1 2]);
                end
                
                lt_plot_zeroline;
                ylim([min(AllSpkDiff_minusDiff_lateminusbase) max(AllSpkDiff_minusDiff_lateminusbase)]);
end



function [MOTIFSTATS_Compiled, SwitchStruct] = lt_neural_v2_ANALY_Swtch_Extract(MOTIFSTATS_Compiled, SwitchStruct, ...
    interpolateCorrNan, RemoveTrialsZeroFR, premotorWind)

%% solution to problem of trials with zero FR (which gives nan for correlation with baselien)

% SOLUTION 1
if interpolateCorrNan==1 & RemoveTrialsZeroFR==0
    % for any trial with fr = 0 througuouit, will get nan for corr with
    % baseline psth. solve by taking avrage of corr of precedinga nd
    % ofllowing trial.
lt_figure; hold on;
lt_plot_text(0, 0.5, 'NOTE: if FR = 0 , neural sim (corr) is interpolated based on flanking trials');
end

% SOLUTION 2 (default)
if RemoveTrialsZeroFR==1
    disp('REMOVING TRIALS WITH 0 FR...');
end

%% lt 7/3/17 - makes sure only one neuron for each chan per switch (no copied data)
% keeps the one with the most data
% also outputs another field will all neurons, including redundant.

%% lt 6/29/17 - for each switch, extract learning stats
% --- plot timecourses and rasters

Numbirds = length(MOTIFSTATS_Compiled.birds);
assert(Numbirds == length(SwitchStruct.bird), 'asdfasdf');

% premotorWind = [-0.075 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset

%%

SwitchStruct.params.premotorWind = premotorWind;

nancounter = 0;
for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
       
       MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
       SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
       
       motiflist = MotifStats.params.motif_regexpr_str;
       nummotifs = length(motiflist);
       exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
       for iii=1:numswitches
          
           assert(length(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron) == ...
               length(MotifStats.neurons), 'asasdf');
           assert(length(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron) == ...
               length(SummaryStruct.birds(1).neurons), 'asasdf');
           
           goodneurons = find([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs] ...
               & [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs]);
           
           swpre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
           swpost = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
           swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
           
           targsyl_tmp = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{1};
           targind = find(strcmp(motiflist, targsyl_tmp));
           
          %% ==== 1) for each motif, plot learning and neural similarity
           for j=1:nummotifs
               
               for nn=goodneurons
                  
                  segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                  
                  
                  if ~isfield(segextract, 'song_filename')
                  continue
                  end
                  
                  %% === remove trials (permanently) with zero FR
                  % within window you care about
                  if RemoveTrialsZeroFR==1
                      tmp = segextract(1).FRsmooth_xbin;
                      premotorInds = find(tmp>(MotifStats.params.motif_predur + premotorWind(1)) ...
                          & tmp<(MotifStats.params.motif_predur + premotorWind(2)));
                      
                      alltrialFR = [segextract.FRsmooth_rate_CommonTrialDur];
                      alltrialFR = alltrialFR(premotorInds, :);
                      
                      trialstoremove = sum(alltrialFR,1)==0;
                      if any(trialstoremove)
                          
                          MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(j).SegmentsExtract(trialstoremove) = [];
                          MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
                          segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                      
                          disp(['REMOVED ' num2str(sum(trialstoremove)) '/' num2str(length(trialstoremove)) '(REASON: fr = 0)']);
                      end
                      
                  end
                  
                  %%
                  % --- extract important things 
                  tmp = segextract(1).FRsmooth_xbin;
                  premotorInds = find(tmp>(MotifStats.params.motif_predur + premotorWind(1)) ...
                    & tmp<(MotifStats.params.motif_predur + premotorWind(2)));
                  
                  MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.premotorInds_FR = premotorInds;
%                 
                  baseInds = [segextract.song_datenum] < swthis & ...
                      [segextract.song_datenum] > swpre;
                  trainInds = [segextract.song_datenum] > swthis & ...
                      [segextract.song_datenum] < swpost;
                  
                  
                  % ============ 1) Neural Similarity (correlation)
                  alltrialFR = [segextract.FRsmooth_rate_CommonTrialDur];
                  alltrialFR = alltrialFR(premotorInds, :);
                  
                  baseFRmean = mean(alltrialFR(:, baseInds),2);
                  
                  NeuralSim = corr(alltrialFR, baseFRmean); % trials x 1
                  
                  if interpolateCorrNan==1
                  if any(sum(alltrialFR,1)==0)
                      % FR =0, will give nan.,. ANY THAT IS NAN WITH THE FLANKING VALUES
                      
                      NeurSimtmp = NeuralSim;
                        if max(find(sum(alltrialFR,1)==0)+1) == length(NeuralSim)+1
                            % then last ind does not have ind+1. replace it
                            % instead with ind-1
                            NeurSimtmp = NeuralSim([1:end, end-1]);
                        end
                            %                       indnext = max([find(sum(alltrialFR,1)==0)+1; length(NeuralSim)])
                      NeuralSim(sum(alltrialFR,1)==0) = ...
                          nanmean([NeurSimtmp(find(sum(alltrialFR,1)==0)-1) NeurSimtmp(find(sum(alltrialFR,1)==0)+1)], 2);
%                           mean([NeuralSim(find(~isnan(NeuralSim))) NeuralSim(find(~isnan(NeuralSim)))],2);
%                       keyboard
                  end
                  end
                  
                  if any(isnan(NeuralSim))
                      disp([birdname '-' exptname '-sw' num2str(iii) '-' motiflist{j}]);
                      disp(NeuralSim)
%                       keyboard
                  end
                  
                  nancounter = nancounter+sum(isnan(NeuralSim));
                  
                  % =========== 2) MEAN FR
                  MeanFR = mean(alltrialFR,1);
                  
                  
                  % =========== 3) STD OF FR (MODULATION)
                  cvFR = std(alltrialFR, 1)./MeanFR;
                  
                  
                  % =====================================================
                  % FF, neural correlation at baseline
                     ffvalstmp = [segextract.FF_val];
                     
                     if sum(baseInds)>2
                     [BaseCorr_FFvsMeanFR, BaseCorr_FFvsMeanFR_pval] = corr(ffvalstmp(baseInds)', MeanFR(baseInds)');
                     [BaseCorr_FFvsNeurSim, BaseCorr_FFvsNeurSim_pval] = corr(ffvalstmp(baseInds)', NeuralSim(baseInds));
                     else
                         BaseCorr_FFvsMeanFR = nan
                         BaseCorr_FFvsMeanFR_pval = nan
                         BaseCorr_FFvsNeurSim = nan
                         BaseCorr_FFvsNeurSim_pval = nan;
                     end
                     
                     
                  % ===========================================================
                  % SNR OF FIRING RATE (baseline, middle, end)
                  
                  
                     
                     
                  % =============================== OUTPUT
                  % time, neural, ff, 
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).NEURvsbase_FRcorr = NeuralSim;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).NEUR_meanFR = MeanFR;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).NEUR_cvFR = cvFR;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).tvals = [segextract.song_datenum];
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).ffvals = [segextract.FF_val];
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds  =baseInds;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds  =trainInds;
                  
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsMeanFR = BaseCorr_FFvsMeanFR;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsMeanFR_pval = BaseCorr_FFvsMeanFR_pval;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsNeurSim = BaseCorr_FFvsNeurSim;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).BASECORR_FFvsNeurSim_pval = BaseCorr_FFvsNeurSim_pval;
                  
                  end
           end
           
           
           % ================== FIGURE OUT WHICH NEURONS ARE REDUNDANT
           channels = [SummaryStruct.birds(1).neurons(goodneurons).channel];
           clusters = [SummaryStruct.birds(1).neurons(goodneurons).clustnum];
           
           numsongs = [];
           neurID = {};
           for j=1:length(goodneurons)
               neurind = goodneurons(j);
               
               tmp = sum(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(neurind).DATA.motif(targind).baseInds) ...
                   + sum(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(neurind).DATA.motif(targind).trainInds);

               numsongs = [numsongs tmp]; % collect number songs
               
               % --- get ID
               neurID = [neurID [num2str(channels(j)) '-' num2str(clusters(j))]];               
           end
           
           % --- if redundant, keep the one with most songs
           uniqueneurons = unique(neurID);
           neuronstoremove = [];
           for j=1:length(uniqueneurons)
               indstmp = find(strcmp(uniqueneurons{j}, neurID));
              if length(indstmp)>1
                  % then multiple version of this dataset. keep the one
                  % with most data
                 
                  [~, maxind] = max(numsongs(indstmp));
                  neuronstoremove = [neuronstoremove indstmp([1:maxind-1 maxind+1:end])];
              end
           end
           
           if ~isempty(neuronstoremove)
               disp('-----');
           disp(['original chans: ' num2str(channels)]);
           disp(['original clusts: ' num2str(clusters)]);
           
           disp(['removed chans: ' num2str([SummaryStruct.birds(1).neurons(goodneurons(neuronstoremove)).channel])]);
           disp(['removed clusts: ' num2str([SummaryStruct.birds(1).neurons(goodneurons(neuronstoremove)).clustnum])]);
           end
           
           goodneurons_all = goodneurons;
           goodneurons(neuronstoremove) = [];
           
           SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons = goodneurons;
           SwitchStruct.bird(i).exptnum(ii).switchlist(iii).goodneurons_redundant = goodneurons_all;
           
       end
    end
end
%     assert(nancounter==0, 'PROBLEM - some corr agaisnt FR of 0, need to interpolate value');
    disp(['--- STILL HAVE ' num2str(nancounter) 'nan corr values in dataset (likely due to no baseline data)']);
    


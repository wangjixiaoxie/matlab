function SwitchStruct = lt_neural_v2_ANALY_Swtch_Extract(MOTIFSTATS_Compiled, SwitchStruct)

%% lt 6/29/17 - for each switch, extract learning stats
% --- plot timecourses and rasters

Numbirds = length(MOTIFSTATS_Compiled.birds);
assert(Numbirds == length(SwitchStruct.bird), 'asdfasdf');

premotorWind = [-0.075 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset

%%

SwitchStruct.params.premotorWind = premotorWind;


for i=1:Numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
       
       MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
       SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
       
       motiflist = MotifStats.params.motif_regexpr_str;
       nummotifs = length(motiflist);
       
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
           
          %% ==== 1) for each motif, plot learning and neural similarity
           for j=1:nummotifs
               
               for nn=goodneurons
                  
                  segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                  
                  if ~isfield(segextract, 'song_filename')
                  continue
                  end
                  
                  % --- extract important things 
                  tmp = segextract(1).FRsmooth_xbin;
                  premotorInds = find(tmp>(MotifStats.params.motif_predur + premotorWind(1)) ...
                    & tmp<(MotifStats.params.motif_predur + premotorWind(2)));
                  
                  baseInds = [segextract.song_datenum] < swthis & ...
                      [segextract.song_datenum] > swpre;
                  trainInds = [segextract.song_datenum] > swthis & ...
                      [segextract.song_datenum] < swpost;
                  
                  
                  % ============ 1) Neural Similarity (correlation)
                  alltrialFR = [segextract.FRsmooth_rate_CommonTrialDur];
                  alltrialFR = alltrialFR(premotorInds, :);
                  
                  baseFRmean = mean(alltrialFR(:, baseInds),2);
                  
                  NeuralSim = corr(alltrialFR, baseFRmean); % trials x 1
                  
                  
                  % =========== 2) MEAN FR
                  MeanFR = mean(alltrialFR,1);
                  
                  
                  % =========== 3) STD OF FR (MODULATION)
                  cvFR = std(alltrialFR, 1)./MeanFR;
                  
                   
                  % =============================== OUTPUT
                  % time, neural, ff, 
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).NEURvsbase_FRcorr = NeuralSim;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).NEUR_meanFR = MeanFR;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).NEUR_cvFR = cvFR;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).tvals = [segextract.song_datenum];
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).ffvals = [segextract.FF_val];
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).baseInds  =baseInds;
                  SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).DATA.motif(j).trainInds  =trainInds;
                  
                  
                 
               end
           end
       end
    end
end



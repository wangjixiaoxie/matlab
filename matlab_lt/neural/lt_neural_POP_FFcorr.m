function MOTIFSTATS_pop = lt_neural_POP_FFcorr(MOTIFSTATS_pop, SummaryStruct)
%% does extraction of data, puts back into struct

%% lt 12/21/17 - firing rates predict pitch? 


%% params

% --- collecting spikes
SpWindow = {'RA', [-0.025 0.025], 'LMAN', [-0.04 0.02]};
SpWindow = {'RA', [-0.015 0.025], 'LMAN', [-0.05 0.02]};
SpWindow = {'RA', [-0.015 0.025], 'LMAN', [-0.05 0.02], 'X', [-0.05 0.02]};
% 
% spwindow = [-0.03 0.03]; % [preonset postonset]

%%

NumBirds = length(MOTIFSTATS_pop.birds);

for i=1:NumBirds
    
    numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);
    params = MOTIFSTATS_pop.birds(i).params;
    
    for ii=1:numexpts
        
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
        for iii=1:numsets
            
            neurons_thisset = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{iii};
            nummotifs = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif);
            DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii);
            nneur = length(neurons_thisset);
            
            for mm = 1:nummotifs
                % ======================= COLLECT MATRIX OF SPIKE COUNTS AND
                % VECTOR OF FF
                if isempty(DAT.motif(mm).SegExtr_neurfakeID)
                    continue
                end
                
                ntrials = length(DAT.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract);
                motifstr = DAT.motif(mm).regexpstr;
                assert(length(DAT.motif(mm).SegExtr_neurfakeID)==length(neurons_thisset), 'asdfas')
                motif_predur = params.motif_predur;
                
                % ------------ 1) FF
                ff = [DAT.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract.FF_val];
                
                if all(isnan(ff))
                    % then this syl does not have defined FF
                    disp(['skipping - FF not defined: ' motifstr]);
                    continue
                end
                
                
                % ------------ 2) Spike count for each neuron
                NspksAll = nan(ntrials, nneur);
                for nn = 1:length(neurons_thisset)
                    % -- which spike window depends on brain region
                    bregion = SummaryStruct.birds(i).neurons(neurons_thisset(nn)).NOTE_Location;
                    if ~any(strcmp({'X', 'LMAN', 'RA'}, bregion))
                        continue
                    end
                    indtmp = find(strcmp(SpWindow, bregion));
                    spwindow = SpWindow{indtmp+1};
                    
                    % --
                    clustnum = SummaryStruct.birds(i).neurons(neurons_thisset(nn)).clustnum;
                    segextract = DAT.motif(mm).SegExtr_neurfakeID(nn).SegmentsExtract;
                    Nspks = lt_neural_QUICK_SpkCnts(segextract, motif_predur, spwindow, clustnum);
                    
                    NspksAll(:, nn) = Nspks;
                end
%                 assert(~any(isnan(NspksAll(:))), 'asdfasd');
                
                % ============================= OUTPUT DATA
                MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST.FF = ff';
                MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).POST.Nspks_allneur = NspksAll;
                                
            end
        end
    end
end


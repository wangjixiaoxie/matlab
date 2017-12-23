function Nspks = lt_neural_QUICK_SpkCnts(segextract, motif_predur, spwindow, clustnum)
%% lt 12/21/17 - extract vector of spike counts, across trials

% segextract is extracted from lt_neural_RegExp (contains data, for a given motif and neuron)
% motif_predur is in sec, onset of locked syl from dat onset
% spwindow = [-0.03 0.03]; % [preonset postonset]
% clustnum optional
%% 

                    ntrials = length(segextract);
                    Nspks = nan(ntrials, 1);
                   
                    if ~exist('clustnum', 'var')
                    assert(sum(unique([segextract.spk_Clust])>0)==1, 'asdfsd'); % i.e. only one cluister...
                    clustnum = [];
                    end
                    windons = motif_predur+spwindow(1);
                    windoff = motif_predur+spwindow(2);
                    for tt = 1:ntrials
                        if isempty(clustnum)
                        indtmp = segextract(tt).spk_Times>windons ...
                            & segextract(tt).spk_Times<windoff;
                            
                        else
                        indtmp = segextract(tt).spk_Clust==clustnum & ...
                            segextract(tt).spk_Times>windons ...
                            & segextract(tt).spk_Times<windoff;
                        end
                        
                        nspks = length(segextract(tt).spk_Times(indtmp));
                        Nspks(tt) = nspks;
                    end

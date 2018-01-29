function OUTSTRUCT = lt_neural_POP_PlotSummary(MOTIFSTATS_pop, SummaryStruct)
%% lt 1/22/18 - extracts data from structure to arrays/cells to prepare for analysis/plotting

plotRaw =1;

%%
NumBirds = length(MOTIFSTATS_pop.birds);

% ============ OUTPUTS
OUTSTRUCT.Pairs.birdID = [];
OUTSTRUCT.Pairs.exptID= [];
OUTSTRUCT.Pairs.motifnum= [];
OUTSTRUCT.Pairs.setnum= [];

OUTSTRUCT.Pairs.neurIDfake = [];
OUTSTRUCT.Pairs.neurIDreal = [];
OUTSTRUCT.Pairs.bregions= {};
OUTSTRUCT.Pairs.bregion_string= {};
OUTSTRUCT.Pairs.ccRealAll= {};
OUTSTRUCT.Pairs.ccShiftAll= {};
OUTSTRUCT.Pairs.ccAuto1_minus= [];
OUTSTRUCT.Pairs.ccAuto2_minus= [];
OUTSTRUCT.Pairs.ccAuto1_mean= [];
OUTSTRUCT.Pairs.ccAuto2_mean= [];
OUTSTRUCT.Pairs.xlags = [];

%% =========== RUN
for i=1:NumBirds
    
    numexpts = length(MOTIFSTATS_pop.birds(i).exptnum);
%     params = MOTIFSTATS_pop.birds(i).params;
    
    for ii=1:numexpts
        
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum);
        for iii=1:numsets
            disp(['bird' num2str(i) '-expt' num2str(ii) '-set' num2str(iii)])
            
            neurons_thisset = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{iii};
            nummotifs = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif);
            %             DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii);
            nneur = length(neurons_thisset);
            
            % =========== WHAT BRAIN REGIONS ARE THESE?
%             BrainregionsList = {SummaryStruct.birds(i).neurons(neurons_thisset).NOTE_Location};
            
            
            
           % ######################### LMAN/RA PAIRS [MODIFY TO BE GENERAL]
            %             indsLMAN = strcmp(BrainregionsList, 'LMAN');
            %             indsRA = strcmp(BrainregionsList, 'RA');
            %             indsX = strcmp(BrainregionsList, 'X');
            %
            %             if onlySetsWithRALMANpaired==1
            %                 if ~any(indsLMAN) & ~any(indsRA)
            %                     continue
            %                 end
            %             end
            %
            for mm = 1:nummotifs
                
%                 segextract_for_trialdur = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                motifstr = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).regexpstr;
                motifpredur = MOTIFSTATS_pop.birds(i).params.motif_predur;
                
                %                 % ========== for plotting raw
                %                 figcount=1;
                %                 subplotrows=5;
                %                 subplotcols=3;
                %                 fignums_alreadyused=[];
                %                 hfigs=[];
                                
                %% 1) EXTRACT CROSSCOVAR FOR ALL PAIRS
                for nn=1:nneur
                    for nnn=nn+1:nneur
                        
                        DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).XCov_neurpair(nn, nnn);
                        
                        if isempty(DAT.bregtmp)
                            continue
                        end
                        
                        % --- brain regions
                        bregtmp = DAT.bregtmp;
                        neurIDfake = [nn nnn];
                        neurIDreal = [neurons_thisset(nn) neurons_thisset(nnn)];
                        
                        % --- xcov
                        ccRealAll = DAT.ccRealAll;
                        ccShiftAll = DAT.ccShiftAll;
                        x = DAT.x;
                        
                        % --- autocovariance
                        ccAuto1 = DAT.ccAuto1;
                        ccAuto1shift = DAT.ccAuto1Shift;
                        ccAuto1_minus = mean(ccAuto1,1) - mean(ccAuto1shift,1);
                        ccAuto1_mean = mean(ccAuto1,1);
                                                
                        ccAuto2 = DAT.ccAuto2;
                        ccAuto2shift = DAT.ccAuto2Shift;
                        ccAuto2_minus = mean(ccAuto2,1) - mean(ccAuto2shift,1);
                        ccAuto2_mean = mean(ccAuto2,1);
                        
                        % ========= collect
                        OUTSTRUCT.Pairs.birdID = [OUTSTRUCT.Pairs.birdID; i];
                        OUTSTRUCT.Pairs.exptID= [OUTSTRUCT.Pairs.exptID; ii];
                        OUTSTRUCT.Pairs.motifnum= [OUTSTRUCT.Pairs.motifnum; mm];
                        OUTSTRUCT.Pairs.setnum= [OUTSTRUCT.Pairs.setnum; iii];
                        OUTSTRUCT.Pairs.xlags = [OUTSTRUCT.Pairs.xlags; x];
                        
                        % ======= flip data to put brain regions in alpha order
                        [~, ind] = sort(bregtmp);
                        if ind(1)>ind(2)
                            % then need to flip
                            
                            bregtmp = fliplr(bregtmp);
                            ccRealAll = fliplr(ccRealAll);
                            ccShiftAll = fliplr(ccShiftAll);
                            neurIDfake = fliplr(neurIDfake);
                            neurIDreal = fliplr(neurIDreal);
                            ccAuto2_minus = fliplr(ccAuto2_minus);
                            ccAuto1_minus = fliplr(ccAuto1_minus);
                            ccAuto1_mean = fliplr(ccAuto1_mean);
                            ccAuto2_mean = fliplr(ccAuto2_mean);
                        end
                        
                        OUTSTRUCT.Pairs.neurIDfake = [OUTSTRUCT.Pairs.neurIDfake; neurIDfake];
                        OUTSTRUCT.Pairs.neurIDreal = [OUTSTRUCT.Pairs.neurIDreal; neurIDreal];
                        OUTSTRUCT.Pairs.bregions= [OUTSTRUCT.Pairs.bregions; bregtmp];
                        OUTSTRUCT.Pairs.ccRealAll=[OUTSTRUCT.Pairs.ccRealAll; ccRealAll];
                        OUTSTRUCT.Pairs.ccShiftAll= [OUTSTRUCT.Pairs.ccShiftAll; ccShiftAll];
                        OUTSTRUCT.Pairs.ccAuto1_minus= [OUTSTRUCT.Pairs.ccAuto1_minus; ccAuto1_minus];
                        OUTSTRUCT.Pairs.ccAuto2_minus= [OUTSTRUCT.Pairs.ccAuto2_minus; ccAuto2_minus];
                        OUTSTRUCT.Pairs.bregion_string = [OUTSTRUCT.Pairs.bregion_string; [bregtmp{1} '-' bregtmp{2}]];
                        
                        OUTSTRUCT.Pairs.ccAuto1_mean= [OUTSTRUCT.Pairs.ccAuto1_mean; ccAuto1_mean];
                        OUTSTRUCT.Pairs.ccAuto2_mean= [OUTSTRUCT.Pairs.ccAuto2_mean; ccAuto2_mean];

                        %% == plot raw?
                        if plotRaw ==1
                            
                            
                        end
                            
                    end
                end
            end
        end
    end
end
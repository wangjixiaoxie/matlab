function MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct)
%% note, made sure to work after lin time warp (

%% lt 12/20/17 - outputs population neurons. makes sure are aligned by trials

%% params

MOTIFSTATS_pop = struct;
MOTIFSTATS_pop.SummaryStruct = SummaryStruct;
%% run
Numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:Numbirds
    
    numexpts = length(SummaryStruct.birds(i).exptnum_pop);

    MOTIFSTATS_pop.birds(i).birdname = SummaryStruct.birds(i).birdname;
    MOTIFSTATS_pop.birds(i).params = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params;
    for ii=1:numexpts
        % ============== for each set of neurons
        numsets = length(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons);
        
        MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons = ...
            SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons;
        
        MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_songfiles = ...
            SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles;
        
        MOTIFSTATS_pop.birds(i).exptnum(ii).exptname = ...
            SummaryStruct.birds(i).exptnum_pop(ii).exptname;
        
        
        for iii=1:numsets
            
            NeurThisSet = SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{iii};
            SongsThisSet = SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles{iii};
            
            motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(NeurThisSet(1)).motif_regexpr_str;
            % sanity checyk
            for mm = 2:length(NeurThisSet)
                assert(all(strcmp(motiflist, MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(NeurThisSet(mm)).motif_regexpr_str)));
            end
            
            % ========================= COLLECT POP DATA FOR ALL MOTIFS
            % go thru each motif. then for each motif go thru each neuron and collect all data that
            % fall within these songfiles
            DATTMP = struct;
            for mm = 1:length(motiflist)
                MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).regexpstr = motiflist{mm};
                
                
                for nn=NeurThisSet
                    
                    % -- skip if no data for this motif
                    if length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif)<mm
                        DATTMP.motif(mm).neur = [];
                        continue
                    end
                    % ---- if this neuron doesn't have data, then skip this
                    % motif
                    if isempty(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract)
                        DATTMP.motif(mm).neur = [];
                        continue
                    end
                    
                    % --- save the rends with correct song filename
                    allfnames = {MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract.song_filename};
                    correctfnames = SongsThisSet;
                    rendstokeep = [];
                    for ff = 1:length(allfnames)
                        thisfname = allfnames{ff};
                        
                        indtmp = findstr(thisfname, '.rhd');
                        thisfname = thisfname(indtmp-13:indtmp-1); % convert to date_time;
                        disp(thisfname);
                        
                        rendstokeep = [rendstokeep any(strcmp(thisfname, correctfnames))];
                    end
                    rendstokeep = logical(rendstokeep);
                    
                    % ===================== EXTRACT DATA TO KEEP (I.E.
                    % WITHIN SHARED SONG FILES)
                    dattokeep = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract(rendstokeep);
                    
                    DATTMP.motif(mm).neur(nn).SegmentsExtract = dattokeep;
                    
                end
            end
            
            % ================================== CHECK THAT ALL DATA
            % ARE ALIGNED ACROSS NEURONS
            for mm = 1:length(motiflist)
                
                if isempty(DATTMP.motif(mm).neur)
                    % then this motif has no data .... SKIP
                    MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID = [];
                    continue
                end
                
                % ---- 1) determine sample size (trials)
                sampsizeall = [];
                for nn=NeurThisSet
                    sampsize = length(DATTMP.motif(mm).neur(nn).SegmentsExtract);
                    sampsizeall = [sampsizeall sampsize];
                end
                assert(length(unique(sampsizeall))==1, 'diff neurons shuld have same N...');
                Ntrials = unique(sampsizeall);
                
                if Ntrials==0
                    % --- then skip this motif
                    MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID = [];
                    continue
                end
                
                % ------- 2) Check that each trial aligned across
                % nuerons AND Collect neuron-specific data.
                for t = 1:Ntrials
                    
                    % -- things I want to compare
                    allmotifdur = [];
                    allfs = [];
                    allsylontimes = [];
                    alldatdur = [];
                    allmsylon = [];
                    allmsyloff = [];
                    for nn=NeurThisSet
                        
                        motifsylonset = DATTMP.motif(mm).neur(nn).SegmentsExtract(t).motifsylOnsets;
                        allmsylon = [allmsylon; motifsylonset];
                        
                        motifsyloffset = DATTMP.motif(mm).neur(nn).SegmentsExtract(t).motifsylOffsets;
                        allmsyloff = [allmsyloff; motifsyloffset];
                        
                        motifdur = DATTMP.motif(mm).neur(nn).SegmentsExtract(t).global_offtime_motifInclFlank ...
                            - DATTMP.motif(mm).neur(nn).SegmentsExtract(t).global_ontime_motifInclFlank;
                        
                        alldatdur = [alldatdur; ...
                            motifdur];
                        
                        allmotifdur = [allmotifdur; ...
                            DATTMP.motif(mm).neur(nn).SegmentsExtract(t).actualmotifdur];
                        
                        allfs= [allfs; ...
                            DATTMP.motif(mm).neur(nn).SegmentsExtract(t).fs];
                        
                        tmp =(DATTMP.motif(mm).neur(nn).SegmentsExtract(t).sylOnTimes_RelDataOnset);
                        tmp = tmp(tmp>0.00001 & tmp<motifdur -0.00001);
                        allsylontimes= [allsylontimes; tmp];
                        
                        %                             DATTMP.motif(mm).neur(nn).SegmentsExtract(t)
                        %
                    end
                    tmp = diff(allmotifdur);
                    assert(all(tmp<0.00001), 'sdasf');
                    tmp = diff(alldatdur);
                    assert(all(tmp<0.00001), 'sdasf');
                    assert(length(unique(allfs))==1, 'sdasf');
                    tmp = diff(allsylontimes,1,1);
                    assert(all(tmp(:)<0.00001), 'sdasf');
                    tmp = diff(allmsylon, 1, 1);
                    assert(all(tmp(:)<0.00001), 'asdfs');
                    tmp = diff(allmsyloff, 1, 1);
                    assert(all(tmp(:)<0.00001), 'asdfs');
                    
                end
                
                % -------------- Take segextract for first neuron (all
                % are identical). Take neural data for all neurons
                %                 segextract_Common = DATTMP.motif(mm).neur(NeurThisSet(1)).SegmentsExtract;
                %                 % --- slide spktimes for all units into this
                %
                %
                %                 for nn=NeurThisSet
                %                         for ttt = 1:Ntrials
                %                             segextract_Common(ttt).spk_Timesnn(nn) = ...
                %                                 DATTMP.motif(mm).neur(nn).SegmentsExtract(ttt).spk_Times;
                %                             segextract(ttt).spk_Clust = ...
                %                                 DATTMP.motif(mm).neur(nn).SegmentsExtract(ttt).spk_Clust;
                %                         end
                %                 end
                
                
                for nn=NeurThisSet
                    segextract = struct;
                    segextract(Ntrials).spk_Times = [];
                    segextract(Ntrials).spk_Clust = [];
                    if nn ==NeurThisSet(1)
                        % if first neuron, keep all info
                        segextract = DATTMP.motif(mm).neur(nn).SegmentsExtract;
                    else
                        
                        for ttt = 1:Ntrials
                            segextract(ttt).spk_Times = ...
                                DATTMP.motif(mm).neur(nn).SegmentsExtract(ttt).spk_Times;
                            segextract(ttt).spk_Clust = ...
                                DATTMP.motif(mm).neur(nn).SegmentsExtract(ttt).spk_Clust;
                        end
                        
                    end
                    
                    %% FINAL OUTPUT
                    indtmp = find(NeurThisSet==nn);
                    MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(indtmp).SegmentsExtract ...
                        = segextract;
                    MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(iii).motif(mm).SegExtr_neurfakeID(indtmp).neurID_orig = nn;
                    
                end
                
                
            end
        end
    end
end



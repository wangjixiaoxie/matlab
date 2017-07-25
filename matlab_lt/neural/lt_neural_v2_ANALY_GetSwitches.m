function [MOTIFSTATS_Compiled, SwitchStruct]= lt_neural_v2_ANALY_GetSwitches(MOTIFSTATS_Compiled)
%% NOTE;
% will take any contingency that is Al (i.e. all WN) and code as same as WN
% off. Reason: to code changes in contingency, should htink iof 100%WN to
% WN (up escapes) and +1.

%%
% If have one neuron that overlaps multiple switches is fine - will extract
% time bins that appropriate for each switch.


%%

LearnSummary = lt_neural_v2_LoadLearnMetadat;



%%
NumBirds = length(MOTIFSTATS_Compiled.birds);
SwitchStruct = struct;

for i=1:NumBirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    birdind = strcmp({LearnSummary.bird.birdname}, birdname);
    learnmetadat = LearnSummary.bird(birdind);
    
%     ExptList = unique([learnmetadat.info(1,:)]);
%     
    ExptList = unique({MOTIFSTATS_Compiled.birds(i).exptnum.exptname});
    
    for ii=1:length(ExptList)
        
        exptname = ExptList{ii};
        
        SummaryStruct = MOTIFSTATS_Compiled.birds(i).exptnum(ii).SummaryStruct;
        assert(length(SummaryStruct.birds)==1, 'asdasdf');
        
        if ~any(strcmp(exptname, {SummaryStruct.birds(1).neurons.exptID}))
            continue
        end
        
        motifstats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        inds = strcmp([learnmetadat.info(1,:)], exptname);
        TargSyls = learnmetadat.info(2,inds);
        
        switchtimes = learnmetadat.info(3:end, inds);
        switchtimes = switchtimes(~cellfun('isempty', switchtimes));
        switchtimes = switchtimes(1:end)'; % become 1d
        %         contingencies = {};
        for j=1:length(switchtimes)
            switchtimes{j} = switchtimes{j}(1:14);
            %             contingencies{j} = switchtimes{j}(16:20);
        end
        disp(switchtimes)
        
        switchtimes = datenum(switchtimes, 'ddmmmyyyy-HHMM');
        
        [~, indtmp] = unique(switchtimes);
        switchtimes = switchtimes(indtmp); % finally, unique switchtimes, over all potentail target syls
        %         contingencies = contingencies(indtmp);
        
        for iii=1:length(switchtimes)
            
            % ================ Figure out which neurons are in time window
            % for this switch
            
            swtime_this = switchtimes(iii);
            if iii==1
                swtime_pre = 0; %
            else
                swtime_pre = switchtimes(iii-1);
            end
            if iii == length(switchtimes)
                swtime_post = switchtimes(iii)+2000;
            else
                swtime_post = switchtimes(iii+1); % 1000 days after ...
            end
            
            neuronlist = find(strcmp({SummaryStruct.birds(1).neurons.exptID}, exptname));
            assert(all(strcmp({SummaryStruct.birds(1).neurons.exptID}, exptname)));
            assert(length(motifstats.neurons)==length(SummaryStruct.birds(1).neurons));
            
            for nn =neuronlist
                
                songdatenums = SummaryStruct.birds(1).neurons(nn).Filedatenum_unsorted;
                %                  = datenum(SummaryStruct.birds(i).neurons(nn).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM');
                
                % ------ ask whether this neuron has files i) before, ii)
                % after, switch (but within bounds of surrounding switches
                
                if any(songdatenums>swtime_pre & songdatenums<swtime_this)
                    haspresongs = 1;
                else
                    haspresongs = 0;
                end
                
                if any(songdatenums>swtime_this & songdatenums<swtime_post)
                    haspostsongs = 1;
                else
                    haspostsongs = 0;
                end
                
                % -- sanity check, if designated to this swich, then
                % shoudl ghave file pre and post.
                if swtime_this == ...
                        datenum(SummaryStruct.birds(1).neurons(nn).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM');
                    assert(haspresongs==1 & haspostsongs==1, 'asfasdfasfasfsdaf');
                end
                
                if (0)   % SANITY check, plot song tdates for this neuron + this switcht ime's borders
                    lt_figure; hold on;
                    title('circles - songs for this neuron; bk: this neuron; bu: learningmetadat');
                    plot(songdatenums, 1, 'ok');
                    plot(swtime_this, 2, 'sb');
                    line([swtime_pre swtime_post], [2 2])
                    line([datenum(SummaryStruct.birds(i).neurons(nn).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM') ...
                        datenum(SummaryStruct.birds(i).neurons(nn).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM')], ylim, ...
                        'Color', 'k')
                    for zzz = 1:length(SummaryStruct.birds(i).neurons(nn).LEARN_WNotherImportantDates)
                        line([datenum(SummaryStruct.birds(i).neurons(nn).LEARN_WNotherImportantDates{zzz}, 'ddmmmyyyy-HHMM') ...
                            datenum(SummaryStruct.birds(i).neurons(nn).LEARN_WNotherImportantDates{zzz}, 'ddmmmyyyy-HHMM')], ...
                            ylim, 'Color', 'k');
                    end
                    ylim([0 3]);
                    xlim([min([songdatenums swtime_this]) max([swtime_this songdatenums])])
                    
                    lt_plot_annotation(1, ['has pre: ' num2str(haspresongs)], 'r');
                    lt_plot_annotation(4, ['has post: ' num2str(haspostsongs)], 'r');
                    pause
                end
                
                % ===== OUTPUT
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).haspresongs = ...
                    haspresongs;
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).haspostsongs = ...
                    haspostsongs;
                
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).exptname = ...
                    SummaryStruct.birds(1).neurons(nn).exptID; % --- neuron in expt?
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).chan = ...
                    SummaryStruct.birds(1).neurons(nn).channel;
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron(nn).dirname = ...
                    SummaryStruct.birds(1).neurons(nn).dirname;
            end
            
            
            % ================ WHAT IS DIRECTION OF LEARNING FOR THE
            % TARGETS IN THIS SWITCH?
            SwtchConts = {}; % {targ, [pre post], ...} where -1, 0, 1 are dn, off, up
            for j = 1:length(TargSyls)
                
                targsyl = TargSyls{j};
                preCont = [];
                postCont = [];
                
                indtmp = strcmp([learnmetadat.info(1,:)], exptname) ... % find this targ syl and this switch time
                    & strcmp([learnmetadat.info(2,:)], targsyl);
                
                % find this transition
                trans_tmp = learnmetadat.info(3:end, indtmp);
                trans_tmp = trans_tmp(~cellfun('isempty', trans_tmp));
                count = 0;
                for jj=1:length(trans_tmp)
                    if swtime_this == datenum(trans_tmp{jj}(1:14), 'ddmmmyyyy-HHMM')
                        % get switch pre and post
                        contingency = trans_tmp{jj}(16:20);
                        
                        if strcmp(contingency(1:2), 'Dn')
                            preCont = -1;
                        elseif strcmp(contingency(1:2), 'Of') | strcmp(contingency(1:2), 'Al')
                            preCont = 0;
                        elseif strcmp(contingency(1:2), 'Up')
                            preCont = 1;
                        end
                        
                        if strcmp(contingency(4:5), 'Dn')
                            postCont = -1;
                        elseif strcmp(contingency(4:5), 'Of') | strcmp(contingency(4:5), 'Al')
                            postCont = 0;
                        elseif strcmp(contingency(4:5), 'Up')
                            postCont = 1;
                        end
                        
                        count = count+1; % sanity check,
                    end
                end
                assert(count ==1, 'should be exactly one for this target')
                assert(~isempty(preCont) & ~isempty(postCont), 'why???');
                
                SwtchConts = [SwtchConts targsyl [preCont postCont]];
            end
            
            % =========== output
            SwitchStruct.bird(i).birdname = birdname;
            SwitchStruct.bird(i).exptnum(ii).exptname = exptname;
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum = switchtimes(iii);
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous = swtime_pre;
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next = swtime_post;
            
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies = SwtchConts;
        end
    end
end

%% ====== extract ff values (for all syls) and ifor about targ

assert(length(MOTIFSTATS_Compiled.birds) == length(SwitchStruct.bird), 'numbirds dont match')

% NOTE: currently only works if only one target syl. Need to modify to
% include if both target syls are trained in same direction

NumBirds = length(SwitchStruct.bird);
for i=1:NumBirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        motifstats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        
        for iii=1:numswitches
            
            
            goodneurons = find([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs] ...
                & [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs]);
            swtime_pre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swtime_post = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swtime_this = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
%             if i==2 & ii==1 && iii==3
%                 keyboard
%             end
            
            % ==== FOR EACH TARG SYL FIGURE OUT DIRECTION OF LEARNING 
            targsyls = motifstats.params.TargSyls;
            learningDirs = {};
            SylLower = {};
            for j=1:length(targsyls)
                
                tsyl = targsyls{j};
                
                indtmp = find(strcmp([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies], tsyl));
                assert(length(indtmp) ==1, 'asdfsdf');
                
                
                if diff(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{indtmp+1})>0;
                    targlearndir = 1;
                elseif diff(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{indtmp+1})<0;
                    targlearndir = -1;
                else
                    targlearndir =0;
                end
                
                learningDirs = [learningDirs tsyl targlearndir];
                
                % ------ syl lower same for both targs?
                syllower = motifstats.params.SingleSyls_inorder(strcmp(motifstats.params.motif_regexpr_str, tsyl));
                SylLower = [SylLower syllower];
            end
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs = learningDirs;
            
            assert(length(unique(SylLower))>0, 'sdasfasdfsdfasd');
            if length(unique(SylLower)) ==1
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl =1;
            else 
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl =0;
            end
               
            
            
            %
            %             if length(targsyls)>1
            %                 SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_targ_ffvals = nan;
            %                 SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_targ_tvals = nan;
            %
            %             else
            %
            %                 tsyl = targsyls{1};
            %
            %                 % -- figure out direction of learning for targ syl
            %                 %             motifstats.params.
            %                 assert(length(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies)==2, 'why? should only be one targ')
            %
            %                 if diff(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{2})>0;
            %                     targlearndir = 1;
            %                 elseif diff(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{2})<0;
            %                     targlearndir = -1;
            %                 else
            %                     targlearndir =0;
            %                 end
            %
            
            
            %                 % ================= COLLECT FF VALS FOR TARG
            %                 % --- go thru all neurons for this switch and collect ff vals
            %                 % (ensure uniqueness by their times)
            %                 %             assert(length(motifstats.neurons) == length(goodneurons), 'hmmm?');
            %                 ffvals_all = [];
            %                 tvals_all = [];
            %
            %                 for j = goodneurons
            %
            %                     indtmp = strcmp(motifstats.params.motif_regexpr_str, tsyl);
            %                     ffvals = [motifstats.neurons(j).motif(indtmp).SegmentsExtract.FF_val];
            %                     tvals = [motifstats.neurons(j).motif(indtmp).SegmentsExtract.song_datenum];
            %
            %                     % --- keep values within switch epoch
            %                     inds_block = tvals>swtime_pre & tvals<swtime_post;
            %                     tvals = tvals(inds_block);
            %                     ffvals = ffvals(inds_block);
            %
            %
            %                     ffvals_all = [ffvals_all ffvals];
            %                     tvals_all = [tvals_all tvals];
            %                 end
            %
            %                 [~, indtmp] = unique(tvals_all); % - get unique values
            %                 tvals_all = tvals_all(indtmp);
            %                 ffvals_all = ffvals_all(indtmp);
            %
            %                 % ========= SAVE VALUES
            %                 SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_targ_ffvals = ffvals_all;
            %                 SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_targ_tvals = tvals_all;
            %             end
            
            
            
            % =============================== EXTRACT FF VALS FOR ALL SYLS
            syllist = motifstats.params.motif_regexpr_str;
            for ss = 1:length(syllist)
                
                sylthis = syllist{ss};
                
                ffvals_all = [];
                tvals_all = [];
                
                for j = goodneurons
                    
                    indtmp = strcmp(motifstats.params.motif_regexpr_str, sylthis);
                    if ~isfield(motifstats.neurons(j).motif(indtmp).SegmentsExtract, 'FF_val')
                    continue
                    end
                    ffvals = [motifstats.neurons(j).motif(indtmp).SegmentsExtract.FF_val];
                    tvals = [motifstats.neurons(j).motif(indtmp).SegmentsExtract.song_datenum];
                    
                    % --- keep values within switch epoch
                    inds_block = tvals>swtime_pre & tvals<swtime_post;
                    tvals = tvals(inds_block);
                    ffvals = ffvals(inds_block);
                    
                    
                    ffvals_all = [ffvals_all ffvals];
                    tvals_all = [tvals_all tvals];
                end
                
                [~, indtmp] = unique(tvals_all); % - get unique values
                tvals_all = tvals_all(indtmp);
                ffvals_all = ffvals_all(indtmp);
                
                % ========= SAVE VALUES
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(ss).ffvals = ffvals_all;
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(ss).tvals = tvals_all;
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(ss).sylname = sylthis;
                
            end
            
            
        end
    end
end

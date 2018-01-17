function [SwitchStruct]= lt_neural_LEARN_getswitch(SummaryStruct)
%% lt 1/16/18 - modified from lt_neural_v2_ANALY_GetSwitches
% NOW DOES NOT NEED MOTIFSTATS


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
NumBirds = length(SummaryStruct.birds);
SwitchStruct = struct;

for i=1:NumBirds

    birdname = SummaryStruct.birds(i).birdname;
    birdind = strcmp({LearnSummary.bird.birdname}, birdname);
    learnmetadat = LearnSummary.bird(birdind);
    ExptList = unique({SummaryStruct.birds(i).neurons.exptID});
    
    for ii=1:length(ExptList)
        exptname = ExptList{ii};

        % ---- find this expt metadat
        inds = strcmp([learnmetadat.info(1,:)], exptname);
        TargSyls = learnmetadat.info(2,inds);
        
        switchtimes = learnmetadat.info(3:end, inds);
        switchtimes = switchtimes(~cellfun('isempty', switchtimes));
        switchtimes = switchtimes(1:end)'; % become 1d

        for j=1:length(switchtimes)
            switchtimes{j} = switchtimes{j}(1:14);
            %             contingencies{j} = switchtimes{j}(16:20);
        end
        disp(switchtimes)
        
        switchtimes = datenum(switchtimes, 'ddmmmyyyy-HHMM');
        
        [~, indtmp] = unique(switchtimes);
        switchtimes = switchtimes(indtmp); % finally, unique switchtimes, over all potentail target syls

        
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
            
            neuronlist = find(strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname));
            
%             neuronlist = find(strcmp({SummaryStruct.birds(1).neurons.exptID}, exptname));
%             assert(all(strcmp({SummaryStruct.birds(1).neurons.exptID}, exptname)));
%             assert(length(motifstats.neurons)==length(SummaryStruct.birds(1).neurons));
            
            for nn=1:length(neuronlist)
                
                neurID = neuronlist(nn);
                
                songdatenums = SummaryStruct.birds(i).neurons(neurID).Filedatenum_unsorted;
                
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
                        datenum(SummaryStruct.birds(i).neurons(neurID).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM')
                    assert(haspresongs==1 & haspostsongs==1, 'asfasdfasfasfsdaf');
                end
                
                if (0)   % SANITY check, plot song tdates for this neuron + this switcht ime's borders
                    lt_figure; hold on;
                    title('circles - songs for this neuron; bk: this neuron; bu: learningmetadat');
                    plot(songdatenums, 1, 'ok');
                    plot(swtime_this, 2, 'sb');
                    
                    line([swtime_pre swtime_post], [2 2])
                    
                    line([datenum(SummaryStruct.birds(i).neurons(neurID).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM') ...
                        datenum(SummaryStruct.birds(i).neurons(neurID).LEARN_WNonDatestr, 'ddmmmyyyy-HHMM')], ylim, ...
                        'Color', 'k')
                    
                    for zzz = 1:length(SummaryStruct.birds(i).neurons(neurID).LEARN_WNotherImportantDates)
                        line([datenum(SummaryStruct.birds(i).neurons(neurID).LEARN_WNotherImportantDates{zzz}, 'ddmmmyyyy-HHMM') ...
                            datenum(SummaryStruct.birds(i).neurons(neurID).LEARN_WNotherImportantDates{zzz}, 'ddmmmyyyy-HHMM')], ...
                            ylim, 'Color', 'k');
                    end
                    ylim([0 3]);
                    xlim([min([songdatenums swtime_this]) max([swtime_this songdatenums])])
                    
                    lt_plot_annotation(1, ['has pre: ' num2str(haspresongs)], 'r');
                    lt_plot_annotation(4, ['has post: ' num2str(haspostsongs)], 'r');
                    pause
                end
                
                % ===== OUTPUT
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neurFakeID(nn).haspresongs = ...
                    haspresongs;
                
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neurFakeID(nn).haspostsongs = ...
                    haspostsongs;
                
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neurFakeID(nn).exptname = ...
                    SummaryStruct.birds(i).neurons(neurID).exptID; % --- neuron in expt?
                
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neurFakeID(nn).chan = ...
                    SummaryStruct.birds(i).neurons(neurID).channel;
                
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neurFakeID(nn).dirname = ...
                    SummaryStruct.birds(i).neurons(neurID).dirname;                
                
                SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neurFakeID(nn).neurID = ...
                    neurID;
                
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


%% ====== extract DIR OF LEARNIGN
NumBirds = length(SwitchStruct.bird);
for i=1:NumBirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);

        for iii=1:numswitches
            
            swtime_pre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swtime_post = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;

            % ==== FOR EACH TARG SYL FIGURE OUT DIRECTION OF LEARNING 
            targsyls = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies(1:2:end);
            
            learningDirs = {};
            SylLower = {};
            for j=1:length(targsyls)
                
                tsyl = targsyls{j};
                
                indtmp = find(strcmp([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies], tsyl));
                assert(length(indtmp) ==1, 'asdfsdf');
                
                if diff(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{indtmp+1})>0
                    targlearndir = 1;
                elseif diff(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningContingencies{indtmp+1})<0
                    targlearndir = -1;
                else
                    targlearndir =0;
                end
                
                learningDirs = [learningDirs tsyl targlearndir];
                
                % ------ syl lower same for both targs?
                syllower = tsyl(strfind(tsyl, '(')+1);
                                
                SylLower = [SylLower syllower];
            end
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs = learningDirs;
            
            assert(length(unique(SylLower))>0, 'sdasfasdfsdfasd');
            if length(unique(SylLower)) ==1
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl =1;
            else 
            SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl =0;
            end
               
            SylLower_targ = SylLower;
            
        end
    end
end

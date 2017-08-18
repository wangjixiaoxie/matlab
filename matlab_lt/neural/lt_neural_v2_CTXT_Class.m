function CLASSIFIEROUT = lt_neural_v2_CTXT_Class(MOTIFSTATS_Compiled, SummaryStruct, nmin)

%% lt 8/11/17

frbinsize = 0.01; % size of bins to reduce Dim of FR vec

%% ===

numbirds = length(MOTIFSTATS_Compiled.birds);

%% === EXTRACT raw FR [don't exclude trials with 0 FR];

RemoveTrialsZeroFR = 0;
% premotorWind = [-0.075 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset
[MOTIFSTATS_Compiled] = lt_neural_v2_ANALY_GetAllFR(MOTIFSTATS_Compiled, ...
    RemoveTrialsZeroFR);

%% === for each neuron/branch point, construct and cross-validate a model

CLASSIFIEROUT = struct;

for i=1:numbirds
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    
    CLASSIFIEROUT.birds(i).birdname = birdname;
    neuroncount = 0;
    for ii=1:numexpts
        
        MotifStats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
        exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        
        numsylctxt = length(MotifStats.params.motif_regexpr_str);
        numneurons = length(MotifStats.neurons);
        %         plotcols = lt_make_plot_colors(numneurons, 0, 0);
        
        
        % ---- go by order of single syls
        SingleSyls = MotifStats.params.singlesyls_unique;
        %         plotcols_single = lt_make_plot_colors(length(SingleSyls), 0, 0);
        %         hsplots = [];
        %         ymax = [];
        
        motif_predur = MotifStats.params.motif_predur;
        
        CLASSIFIEROUT.birds(i).exptnum(ii).exptname = exptname;
        CLASSIFIEROUT.birds(i).exptnum(ii).sylctxtlist = MotifStats.params.motif_regexpr_str;
        
        
        for ss = 1:length(SingleSyls)
            singlesyl = SingleSyls{ss};
            
            CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).singlesyl = singlesyl;
            
            Xall = []; % trials x bins (FR vectors)
            xtimesall = []; % 1 x bins (time stamps)
            Neuron = []; % trials x 1
            Y = []; % context indicator
            
            CtxtClasses = {};
            contextcounter = 0;
            % ---
            for j=1:numsylctxt
                
                sylname = MotifStats.params.motif_regexpr_str{j};
                
                % -- only continue if syl corresponds to this single syl
                syltmp = sylname(strfind(sylname, '(')+1);
                if ~strcmp(singlesyl, syltmp)
                    continue
                end
                
                contextcounter = contextcounter+1;
                CtxtClasses = [CtxtClasses sylname];
                for nn=1:numneurons
                    
                    segextract = MotifStats.neurons(nn).motif(j).SegmentsExtract;
                    clustnum = MotifStats.neurons(nn).clustnum;
                    %                     chan = MotifStats.neurons(nn).motif(j).Params.channel_board;
                    
                    if ~isfield(segextract, 'FRsmooth_rate')
                        % then no data
                        continue
                    end
                    
                    
                    % ==
                    
                    % 1) ----------- figure out inds for this time window
                    frtimewindow = [-0.075 0.025]; % on and off, relative to syl onset
                    xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    indsFR = xbin>=(motif_predur+frtimewindow(1)) & xbin<=(motif_predur+frtimewindow(2));
                    
                    X = [segextract.FRsmooth_rate_CommonTrialDur];
                    X = X(indsFR, :);
                    X = X';
                    
                    xtimes = xbin(indsFR);
                    xtimes = xtimes';
                    
                    % --------------- reduce Dim of vectors (by binning)
                    TrimDown = 1;
                    [X, xtimes] = lt_neural_v2_QUICK_binFR(X, xtimes, frbinsize, TrimDown);
                    
                    Xall = [Xall; X];
                    xtimesall = [xtimesall; xtimes];
                    
                    % ----- output context and neurons
                    Neuron = [Neuron; nn*ones(size(X,1),1)];
                    Y = [Y; contextcounter*ones(size(X,1),1)];
                    
                    
                    %                     % --- PLOT SMOOTHED FR
                    %                     xbin = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    %                     motif_predur = MotifStats.params.motif_predur;
                    %                     inds = xbin>(motif_predur+WindowToPlot(1)) & xbin<(motif_predur+WindowToPlot(2));
                    %
                    %                     FRall = [segextract.FRsmooth_rate_CommonTrialDur];
                    %
                    
                end
            end
            
            % --- convert Y to categorical array (with names of contexts)
            Y = categorical(Y);
            
            %% ================= MODEL
            
            % =================== 1 model for each neuron
            neuronIDs = unique(Neuron);
            didgetdata = 0;
            for nn = neuronIDs'
                
                Xthisneuron = Xall(Neuron ==nn, :);
                Ythisneuron = Y(Neuron ==nn);
                
                ctxtlist = unique(Ythisneuron);
                
                if length(ctxtlist)<2
                    % then only one context, skip this
                    continue
                end
                
                
                
                % TO DO: make sure each context has at least N datapoints
                
                Nall = grpstats(1:length(Ythisneuron), Ythisneuron, 'numel');
                
                % - for any context that has <N datapts, remove it
                if any(Nall < nmin)
                    indstoremove = [];
                    % figure out which ctxt
                    for zz = ctxtlist'
                        if sum(Ythisneuron==zz) < nmin
                            % remve thiese inds
                            indstoremove = [indstoremove find(Ythisneuron==zz)];
                        end
                    end
                    
                    Xthisneuron(indstoremove, :) = [];
                    Ythisneuron(indstoremove) = [];
                    
                    ctxtlist = unique(Ythisneuron);
                    
                    % -- decide if have to skip now
                    if length(ctxtlist)<2
                        % then only one context, skip this
                        continue
                    end
                end
                
                
                % ====================================== LOO prediction
                [Ypredicted, ConfMat, accuracy, sensitivity_mean] ...
                    = lt_neural_v2_QUICK_classify(Xthisneuron, Ythisneuron);
                
                
                % ================= OUTPUT
                neuroncount = neuroncount+1;
                CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).neuron(nn).ctxts_that_exist = unique(Ythisneuron);
                CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).neuron(nn).ctxts_that_exist_num = length(unique(Ythisneuron));
                CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).neuron(nn).ctxts_that_exist_name = CtxtClasses(unique(Ythisneuron));
                
                CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).neuron(nn).Rslts_ConfMat = ConfMat;
                CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).neuron(nn).Rslts_accuracy = accuracy;
                CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).neuron(nn).Rslts_sensitivity_mean = sensitivity_mean;
                
                disp(['done classify for ' birdname '-' exptname '-' singlesyl '-n' num2str(nn)]);
                
                didgetdata = 1;
            end
            
            if didgetdata==1
               % only save this metadata if have actual data.
               CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).X_frvec = Xall;
            CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).xtimesall_frtbins = xtimesall;
            CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).Neuron = Neuron;
            CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).Y_ctxtID = Y;
            CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(ss).ctxtID_classnames = CtxtClasses;
            
            
            end
            
        end
    end
end




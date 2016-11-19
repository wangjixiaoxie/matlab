function lt_neural_MultNeur_AllOnebackNeur(NeuronDatabase, TRANSMATRIX)
%% by default ignores transitions containing "-"


NumNeurons=length(NeuronDatabase.neurons);

motif_predur=0.13;
motif_postdur=0.05;


spktimefield='spk_Times';
window = 0.02;
windshift = 0.004;

mindatFreq=8; % only analyze transitions with N or more cases
%%

assert(length(TRANSMATRIX.neuron) == length(NeuronDatabase.neurons), 'asdfacadac');



%% ================ CONVERGENT
branchtype = 'conv';
for i=1:NumNeurons
    lt_figure; hold on;
    % ========================= FIRST, extract data
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    
    % =========================== SECOND, for each potential one-back
    % transition, extract motif data and smoothed firing rate
    %     plotcols=lt_make_plot_colors(size(TRANSMATRIX.neuron(i).matrix,1), 0, 0);
    for j=1:length(TRANSMATRIX.neuron(i).sysInOrder)
        syl1=TRANSMATRIX.neuron(i).sysInOrder(j);
        
        if strcmp(branchtype, 'conv')
            plotcol=[rand rand rand];
        end
        for jj=1:length(TRANSMATRIX.neuron(i).sysInOrder)
            syl2=TRANSMATRIX.neuron(i).sysInOrder(jj);
            if strcmp(branchtype, 'div')
                plotcol=[rand rand rand];
            end
            
            % 1) -----------  exceptions
            if TRANSMATRIX.neuron(i).matrix(j,jj) < mindatFreq
                continue
            elseif strcmp(syl1, '-') || strcmp(syl2, '-');
                continue
            end
            
            % 2) ---------------------------------- extract motif data
            if strcmp(branchtype, 'conv')
                regexpr_str=[syl1 '(' syl2 ')']; % token is 2nd syl, locks to that.
            elseif strcmp(branchtype, 'div')
                regexpr_str=['(' syl1 ')' syl2]; %
            end
            
            % --- GET SMOOTHED FR
            disp(regexpr_str)
            
            predur=motif_predur; % sec
            postdur=motif_postdur; % sec
            alignByOnset=1;
            if strcmp(branchtype, 'div')
                alignByOnset=0;
            end
            
            WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
            % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
            
            % 3) ---------------------- get smoothed firing rate over all trials
            numtrials=length(SegmentsExtract);
            clustnum=NeuronDatabase.neurons(i).clustnum;
            Yspks={};
            for k=1:numtrials
                inds=SegmentsExtract(k).spk_Clust==clustnum;
                spktimes=SegmentsExtract(k).(spktimefield)(inds);
                Yspks{k}=spktimes;
            end
            
            % -- convert to smoothed rate
            [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            % --- GET SMOOTHED FR [end]

            
            % 4) --------- PLOT in correct position
            if strcmp(branchtype, 'conv')
                lt_subplot(8, 3, jj); hold on; % convergent
                title(['to ' syl2]);
            elseif strcmp(branchtype, 'div')
                lt_subplot(8, 3, j); hold on; %
                title(['from ' syl1]);
            end
            
            shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
            if strcmp(branchtype, 'conv')
                lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl1, plotcol);
            elseif strcmp(branchtype, 'div')
                lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl2, plotcol);
            end
            
            % 6) -------- final things
            line([predur predur], ylim, 'Color', 'k');
            axis tight
            lt_plot_zeroline;
            
            
            % 5) ------- overlay timing of each syl (the token syl)
            if strcmp(branchtype, 'conv')
                
                token_offtime = median(SongDat.AllOffsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                line([predur token_offtime], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
                
            elseif strcmp(branchtype, 'div')
                token_ontime = median(SongDat.AllOnsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                line([token_ontime predur], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
            end
            
            tmp=get(gca, 'Ylim');
            set(gca, 'Ylim', [-16 tmp(2)]);
            
        end
    end
    if strcmp(branchtype, 'conv')
        lt_subtitle(['CONV - neuron ' num2str(i)]);
    elseif strcmp(branchtype, 'div')
        lt_subtitle(['DIV - neuron ' num2str(i)]);
    end
end


%% ===== SAME, BUT DIVERGENT
branchtype = 'div';
for i=1:NumNeurons
    lt_figure; hold on;
    % ========================= FIRST, extract data
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    
    % =========================== SECOND, for each potential one-back
    % transition, extract motif data and smoothed firing rate
    %     plotcols=lt_make_plot_colors(size(TRANSMATRIX.neuron(i).matrix,1), 0, 0);
    for j=1:length(TRANSMATRIX.neuron(i).sysInOrder)
        syl1=TRANSMATRIX.neuron(i).sysInOrder(j);
        
        if strcmp(branchtype, 'conv')
            plotcol=[rand rand rand];
        end
        for jj=1:length(TRANSMATRIX.neuron(i).sysInOrder)
            syl2=TRANSMATRIX.neuron(i).sysInOrder(jj);
            if strcmp(branchtype, 'div')
                plotcol=[rand rand rand];
            end
            
            % 1) -----------  exceptions
            if TRANSMATRIX.neuron(i).matrix(j,jj) < mindatFreq
                continue
            elseif strcmp(syl1, '-') || strcmp(syl2, '-');
                continue
            end
            
            % 2) ---------------------------------- extract motif data
            if strcmp(branchtype, 'conv')
                regexpr_str=[syl1 '(' syl2 ')']; % token is 2nd syl, locks to that.
            elseif strcmp(branchtype, 'div')
                regexpr_str=['(' syl1 ')' syl2]; %
            end
            
            disp(regexpr_str)
            
            predur=motif_predur; % sec
            postdur=motif_postdur; % sec
            alignByOnset=1;
            if strcmp(branchtype, 'div')
                alignByOnset=0;
            end
            
            WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
            % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur);
            
            % 3) ---------------------- get smoothed firing rate over all trials
            numtrials=length(SegmentsExtract);
            clustnum=NeuronDatabase.neurons(i).clustnum;
            Yspks={};
            for k=1:numtrials
                inds=SegmentsExtract(k).spk_Clust==clustnum;
                spktimes=SegmentsExtract(k).(spktimefield)(inds);
                Yspks{k}=spktimes;
            end
            
            % -- convert to smoothed rate
            [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            
            % 4) --------- PLOT in correct position
            if strcmp(branchtype, 'conv')
                lt_subplot(8, 3, jj); hold on; % convergent
                title(['to ' syl2]);
            elseif strcmp(branchtype, 'div')
                lt_subplot(8, 3, j); hold on; %
                title(['from ' syl1]);
            end
            
            shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
            if strcmp(branchtype, 'conv')
                lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl1, plotcol);
            elseif strcmp(branchtype, 'div')
                lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl2, plotcol);
            end
            
            % 6) -------- final things
            line([predur predur], ylim, 'Color', 'k');
            axis tight
            lt_plot_zeroline;
            
            
            % 5) ------- overlay timing of each syl (the token syl)
            if strcmp(branchtype, 'conv')
                
                token_offtime = median(SongDat.AllOffsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                line([predur token_offtime], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
                
            elseif strcmp(branchtype, 'div')
                token_ontime = median(SongDat.AllOnsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                line([token_ontime predur], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
            end
            
            tmp=get(gca, 'Ylim');
            set(gca, 'Ylim', [-16 tmp(2)]);
            
        end
    end
    if strcmp(branchtype, 'conv')
        lt_subtitle(['CONV - neuron ' num2str(i)]);
    elseif strcmp(branchtype, 'div')
        lt_subtitle(['DIV - neuron ' num2str(i)]);
    end
end


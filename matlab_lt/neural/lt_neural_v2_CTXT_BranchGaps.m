function ALLBRANCH = lt_neural_v2_CTXT_BranchGaps(ALLBRANCH)
%% lt 9/27/17 - get anovas on syl and gap durs (i.e. across classes for a given branch point


%% for each branch, calculate an anova on gap and syl durations

numalign = length(ALLBRANCH.alignpos);

AllSylOmega = [];
AllGappreOmega = [];
AllGappostOmega = [];


for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    motifpredur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpredur;
    
    for ii=1:numbirds
        
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        numbranch = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        for j=1:numbranch
            numneuron = length(ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron);
            
            for nn=1:numneuron
                
                datneur = ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn);
                
                if isempty(datneur.SylGapDurs)
                    continue
                end
                
                % ############################################# SYL
                X = {datneur.SylGapDurs.classnum.Dur_syl};
                Xout = [];
                Levels = [];
                for k=1:length(X)
                    Xout = [Xout X{k}];
                    Levels = [Levels; k*ones(length(X{k}),1)];
                end
                
                [p, tabletmp] = anova1(Xout, Levels, 'off');
                
                ss_effect = tabletmp{2,2};
                ss_total = tabletmp{4,2};
                df_effect = tabletmp{2,3};
                ms_error = tabletmp{3,4};
                
                
                % calculate omega squared
                numerator = ss_effect - df_effect*ms_error;
                denominator = ss_total + ms_error;
                omega2 = numerator/denominator;
                
                
                % ================================= PLUG B ACK IN
                ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn).DurAnovas.syl_omega = omega2;
                ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn).DurAnovas.syl_p = p;
                AllSylOmega = [AllSylOmega omega2];

                
                % ############################################# PREGAP
                X = {datneur.SylGapDurs.classnum.Dur_gappre};
                Xout = [];
                Levels = [];
                for k=1:length(X)
                    Xout = [Xout X{k}];
                    Levels = [Levels; k*ones(length(X{k}),1)];
                end
                
                [p, tabletmp] = anova1(Xout, Levels, 'off');
                
                ss_effect = tabletmp{2,2};
                ss_total = tabletmp{4,2};
                df_effect = tabletmp{2,3};
                ms_error = tabletmp{3,4};
                
                
                % calculate omega squared
                numerator = ss_effect - df_effect*ms_error;
                denominator = ss_total + ms_error;
                omega2 = numerator/denominator;
                
                
                % ================================= PLUG B ACK IN
                ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn).DurAnovas.gappre_omega = omega2;
                ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn).DurAnovas.gappre_p = p;
                 AllGappreOmega = [AllGappreOmega omega2];
               
                
                % ############################################# POSTGAP
                X = {datneur.SylGapDurs.classnum.Dur_gappost};
                Xout = [];
                Levels = [];
                for k=1:length(X)
                    Xout = [Xout X{k}];
                    Levels = [Levels; k*ones(length(X{k}),1)];
                end
                
                [p, tabletmp] = anova1(Xout, Levels, 'off');
                
                ss_effect = tabletmp{2,2};
                ss_total = tabletmp{4,2};
                df_effect = tabletmp{2,3};
                ms_error = tabletmp{3,4};
                
                
                % calculate omega squared
                numerator = ss_effect - df_effect*ms_error;
                denominator = ss_total + ms_error;
                omega2 = numerator/denominator;
                
                
                % ================================= PLUG B ACK IN
                ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn).DurAnovas.gappost_omega = omega2;
                ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn).DurAnovas.gappost_p = p;
                 AllGappostOmega = [AllGappostOmega omega2];
            end
            
            
        end
    end
end


%% ==== global distribution of omegas

lt_figure; hold on;

% syl
lt_subplot(2,2,1); hold on;
title('syl');
lt_plot_histogram(AllSylOmega);
xlabel('omega2');

% gappre
lt_subplot(2,2,2); hold on;
title('gap pre');
lt_plot_histogram(AllGappreOmega);
xlabel('omega2');

% syl
lt_subplot(2,2,3); hold on;
title('gap post');
lt_plot_histogram(AllGappostOmega);
xlabel('omega2');



%% in progress, pick an aligntment positiona nd bird, will plot all syl, gap(pre an dpost) durs
i = 1;
ii=1;

lt_figure; hold on;
ylabel('dur (syl:black; pre:green; post:red');
xlabel('branch (neurons for a given branch, overlayed');

numbranch = length(ALLBRANCH.alignpos(i).bird(ii).branch);
count = 0;
for iii=1:numbranch
    numneurons = length(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron);
    
    SylMeans = [];
    for nn=1:numneurons
        
        if isempty(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs)
            continue
        end
        
        if isempty(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn))
            continue
        end
        
        % =============== get means
        functmp = @(x)mean(x);
        SylMeans = cellfun(functmp, ...
            {ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs.classnum.Dur_syl});
        GapPreMeans = cellfun(functmp, ...
            {ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs.classnum.Dur_gappre});
        GapPostMeans = cellfun(functmp, ...
            {ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs.classnum.Dur_gappost});
        
        
        % ============== get std
        functmp = @(x)std(x);
        SylStd = cellfun(functmp, ...
            {ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs.classnum.Dur_syl});
        GapPreStd = cellfun(functmp, ...
            {ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs.classnum.Dur_gappre});
        GapPostStd = cellfun(functmp, ...
            {ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).SylGapDurs.classnum.Dur_gappost});
        
        
        % ======== anova results
        
        % =========== plot summary of all branches for this bird
        x = (1:length(SylMeans)) + count;
        %-- syl
        lt_plot(x, SylMeans, {'Errors', SylStd, 'LineStyle', '-', 'LineWidth', 1}); % syl
        if ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).DurAnovas.syl_p<0.05
            omega = ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).DurAnovas.syl_omega;
            lt_plot_text(x(end), SylMeans(end), ['R2=' num2str(omega)], 'k');
        end
        
        % -- gap pre
        lt_plot(x, GapPreMeans, {'Errors', GapPreStd, 'LineStyle', '-', 'Color', 'g'});
        if ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).DurAnovas.gappre_p<0.05
            omega = ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).DurAnovas.gappre_omega;
            lt_plot_text(x(end), GapPreMeans(end), ['R2=' num2str(omega)], 'g');
        end
        
        
        % -- gap post
        lt_plot(x, GapPostMeans, {'Errors', GapPostStd, 'LineStyle', '-', 'Color', 'r'});
        if ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).DurAnovas.gappost_p<0.05
            omega = ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).DurAnovas.gappost_omega;
            lt_plot_text(x(end), GapPostMeans(end), ['R2=' num2str(omega)], 'r');
        end
        
        %         disp(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nn).prms_regexpstrlist)
        
    end
    count = count+length(SylMeans);
    
end
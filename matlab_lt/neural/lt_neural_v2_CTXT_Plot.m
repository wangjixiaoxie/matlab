function lt_neural_v2_CTXT_Plot(CLASSIFIEROUT)



%% COLLECT classifier accuracy

numbirds = length(CLASSIFIEROUT.birds);

AllNumCtxts = [];
AllAccuracy = [];
AllSensitivity = [];


for i=1:numbirds
   numexpts = length(CLASSIFIEROUT.birds(i).exptnum);
   
   for ii=1:numexpts
      
       numbranches = length(CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum);
       
       for iii=1:numbranches
           
          if ~isfield(CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(iii), 'neuron')
              % then does hnot have data (e.g. not enough contexts, not
              % enough data)
              continue
          end
          
          if isempty(CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(iii).neuron)
              continue
          end
          
          
          % ---------------- COLLECT
          numneurons = length(CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(iii).neuron);
          
          for nn =1:numneurons 
             numctxts = CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(iii).neuron(nn).ctxts_that_exist_num;
             
             accuracy = CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(iii).neuron(nn).Rslts_accuracy;
             sensitivity = CLASSIFIEROUT.birds(i).exptnum(ii).singlesylnum(iii).neuron(nn).Rslts_sensitivity_mean;
             
             % =========== COLLECT
              AllNumCtxts = [AllNumCtxts; numctxts];
 AllAccuracy = [AllAccuracy accuracy];
AllSensitivity = [AllSensitivity sensitivity];
             
          end
             
       
       
       end
           
       
   end
   
end


%% ========= PLOT ALL

lt_figure; hold on;

lt_subplot(2,2,1); hold on;
xlabel('numctxts');
ylabel('accuracy');

plot(AllNumCtxts, AllAccuracy, 'ok');

[ymean, yN] = grpstats(AllAccuracy, AllNumCtxts, {'mean', 'numel'});

ctxtslist = unique(AllNumCtxts);

lt_plot_bar(ctxtslist, ymean, {'Color', 'b'});
lt_plot_bar(ctxtslist, 1./ctxtslist, {'Color', 'k'}); % null

xlim([0 max(AllNumCtxts)+1]);
ylim([0 1]);



lt_subplot(2,2,2); hold on;
xlabel('numctxts');
ylabel('sensitivity (pr(pred A if A))');

plot(AllNumCtxts, AllSensitivity, 'ok');

[ymean, yN] = grpstats(AllSensitivity, AllNumCtxts, {'mean', 'numel'});

ctxtslist = unique(AllNumCtxts);

lt_plot_bar(ctxtslist, ymean, {'Color', 'b'});
lt_plot_bar(ctxtslist, 1./ctxtslist, {'Color', 'k'}); % null

xlim([0 max(AllNumCtxts)+1]);
ylim([0 1]);



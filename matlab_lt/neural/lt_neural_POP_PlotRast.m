
plotCols = lt_make_plot_colors(3, 0, 0);
plotColsRegions = {'LMAN', 'RA', 'X'};

bb =1;
ee =2;
ss = 2;
mm = 1;

% ========= fopr this set
nthis = MOTIFSTATS_pop.birds(bb).exptnum(ee).Sets_neurons{ss};
locthis = {SummaryStruct.birds(bb).neurons(nthis).NOTE_Location};


% ========= for this motif
lt_figure; hold on;
for k=1:length(nthis)
    clustnum = SummaryStruct.birds(bb).neurons(nthis(k)).clustnum;
    dat = MOTIFSTATS_pop.birds(bb).exptnum(ee).DAT.setnum(ss).motif(mm).SegExtr_neurfakeID(k).SegmentsExtract;
%     shiftpos = -0.3 + (k-1)*0.6/(length(nthis)-1);
    loc = locthis(k);
    
    plotcol = plotCols{strcmp(plotColsRegions, loc)};
    
    % --- for each trial, plot
    ntrial = length(dat);
    for t =1:ntrial
       stimes = dat(t).spk_Times(dat(t).spk_Clust == clustnum);
       
%        ypos = t+shiftpos;
       ypos = (t-1)*(length(nthis)+1) + k;
       lt_neural_PLOT_rasterline(stimes, ypos, plotcol, 1);
       
    end
end


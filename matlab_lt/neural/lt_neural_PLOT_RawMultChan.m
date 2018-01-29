function lt_neural_PLOT_RawMultChan(dirname, batchf, chanstoplot);

%% lt 11/5/17 - given a batch file and channel, plots neural voltage trace 
% can overlay multiple channels
% can plot individual songs or entire thing

% dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/NOTSINGING';
% batchf = 'BatchTest';
% chanstoplot = [9 14 17 18 21];

%%
cd(dirname);

    fid = fopen(batchf);
    filename = fgetl(fid);
    while ischar(filename)
        
        
        ChansToPlot.DigChans_zero=[0]; % make string "all" to plot all that exist. empty array to ignore
        ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
        % ChansToPlot.AmpChans_zero=[9 14 19];
        ChansToPlot.AmpChans_zero=chanstoplot;
        % ChansToPlot.AmpChans_zero=[8 9 11 16 17 20 21];
        % ChansToPlot.AmpChans_zero=[8 9 16 17 20 21];
        
        % neuralFiltLow=500;
        neuralFiltLow=300;
        
        PlotWhat.raw=0;
        PlotWhat.filt=1;
        PlotWhat.rect_sm=0;
        PlotWhat.raster=0;
        PlotWhat.digital=0;
        
        Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
        Raster.ThrXNoise=3; % threshold for units, used for all channels, uses absolute data for peak detection
        Raster.PosOrNeg=-1; % -1 or +1, for crossings.
        
        numsubplots = min([length(chanstoplot), 6]);
        if numsubplots==length(chanstoplot)-1
            numsubplots = numsubplots-1;
        end
        
        lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster, numsubplots)
        
        disp([' ======== ' filename]);
        
        pause;
        filename = fgetl(fid);
    end

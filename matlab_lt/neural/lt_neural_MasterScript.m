
%% ===== plot song

lt_neural_extractaudio


%% ====== plot single file dat [align neural and song]

filename='bk7_audioOnly_160807_071018.rhd';
ChansToPlot.DigChans_zero=[]; % make string "all" to plot all that exist. empty array to ignore
ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
ChansToPlot.AmpChans_zero=[9 14 19 23];
neuralFiltLow=500;

PlotWhat.raw=0;
PlotWhat.filt=1;
PlotWhat.rect_sm=1;
PlotWhat.raster=1;

Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.

lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm)
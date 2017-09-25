function [SongDat, NeurDat, Params] = lt_neural_RASamMel_ExtractDat(SummaryStruct, birdnum, neurnum)
%% lt 9/25/17 - extracts SongDat and NeurDat from SummaryStruct

datstruct = SummaryStruct.birds(birdnum).neurons(neurnum).dat;
birdname = SummaryStruct.birds(birdnum).birdname;

%%

SongDat = struct;
SongDat.AllLabels = datstruct.labels;
SongDat.AllOnsets = datstruct.onsets./1000; % ms
SongDat.AllOffsets = datstruct.offsets./1000; % ms
SongDat.FFvals = datstruct.peak_pinterp_labelvec;

NeurDat = struct;
NeurDat.metaDat(1).fs = datstruct.Fs;
NeurDat.TotalSamps = datstruct.length_song;
NeurDat.spiketimes = datstruct.spiketimes;
NeurDat.filenames = datstruct.fname_arr;

Params.birdname = birdname;
Params.exptname = '';



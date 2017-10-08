%% lt 9/24/17 - for analysis of RA data. Will converge with other analyses I wrote

%%  extract summary struct (includes spk data as well)

%% =================== EXTRACT SOBER/WOHL RA DATA
clear all; close all;
SummaryStruct = lt_neural_RASamMel_SummaryStruct;

%% ====================== LIST LABELS AND DATES FOR ALL NEURONS

for i=1:length(SummaryStruct.birds)
    disp(' ');
    disp(' ############################# ');
    disp(SummaryStruct.birds(i).birdname);
    AllSyls = [];
    for ii=1:length(SummaryStruct.birds(i).neurons)
        date = SummaryStruct.birds(i).neurons(ii).date;
        disp([date ':: ' SummaryStruct.birds(i).neurons(ii).dat.labels])
        disp(strfind(SummaryStruct.birds(i).neurons(ii).dat.labels, 'O'))
        AllSyls = [AllSyls SummaryStruct.birds(i).neurons(ii).dat.labels];
    end
    disp(['UNIQUE: ' unique(AllSyls)]);
end

%% ==================== FOR EACH NEURON, PLOT
% 1) SONGS, LABELS; 2) RAW DAT; 3) RASTERS
close all
lt_neural_RASamMel_PlotRaw(SummaryStruct);

%% CONTEXT ANALYSIS - USE CONTEXT IN lt_neural_MasterScript_v2

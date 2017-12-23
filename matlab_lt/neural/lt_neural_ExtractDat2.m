function [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, birdnum, ...
    neurnum, extractSound)
%% note: autmoatically only keeps spikes of this cluster

%%
if ~exist('extractSound', 'var')
    extractSound = 0;
end

if isempty(extractSound)
    extractSound=0;
end

%%

i = birdnum;
ii = neurnum;

if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
    % then is sam/mel ra data
    
    [SongDat, NeurDat, Params] = lt_neural_RASamMel_ExtractDat(SummaryStruct, i, ii);
    
    
else
    
    cd(SummaryStruct.birds(i).neurons(ii).dirname);
    cd ..
    batchf=SummaryStruct.birds(i).neurons(ii).batchfilename;
    channel_board=SummaryStruct.birds(i).neurons(ii).channel;
    clustToKeep =SummaryStruct.birds(i).neurons(ii).clustnum;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, ...
        extractSound, clustToKeep);
    Params.birdname = SummaryStruct.birds(i).birdname;
    Params.exptname = SummaryStruct.birds(i).neurons(ii).exptID;
    Params.neurnum = neurnum;
end
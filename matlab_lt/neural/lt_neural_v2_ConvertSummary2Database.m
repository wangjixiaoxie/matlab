function [NeuronDatabase, SummaryStruct] = ...
    lt_neural_v2_ConvertSummary2Database(BirdsToKeep)

%% TO DO:
% 1) Extract "notes" from summary struct


%% lt 3/6/17 - converts from summary structure to neuron database to use in analyses
% in conversion can filter to only extract certain experiments/birds, etc

% BirdsToKeep = {'bk7'}; % optional - if empty, then gets all birds.


%%
SummaryStruct =  lt_neural_v2_LoadSummary;

if isempty(BirdsToKeep)
    BirdsToKeep = unique({SummaryStruct.birds.birdname}); % keep all birds
end

%% =====

numbirds = length(SummaryStruct.birds);
NeuronDatabase = struct;
count = 1;

for i=1:numbirds
    
    if ~any(strcmp(SummaryStruct.birds(i).birdname, BirdsToKeep))
        disp(['skipped ' SummaryStruct.birds(i).birdname]);
        continue
    end
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    birdname = SummaryStruct.birds(i).birdname;
    
    for ii=1:numneurons
        
        tmpstruct = SummaryStruct.birds(i).neurons(ii);
        disp('---- EXTRACTING TO NEURON DATABASE ...');
        disp(tmpstruct);
        
        NeuronDatabase.neurons(count).birdname = birdname; %
        NeuronDatabase.neurons(count).exptID=tmpstruct.exptID; %
        NeuronDatabase.neurons(count).date=tmpstruct.date; % date
        NeuronDatabase.neurons(count).batchfile=tmpstruct.batchfilename; % batchfile (songs)
        NeuronDatabase.neurons(count).chan=tmpstruct.channel; % channel
        NeuronDatabase.neurons(count).clustnum=tmpstruct.clustnum; % cluster
        NeuronDatabase.neurons(count).electrode_depth=tmpstruct.electrode_depth; % cluster
        NeuronDatabase.neurons(count).NOTE_is_single_unit=tmpstruct.NOTE_is_single_unit; % is_single_unit
        NeuronDatabase.neurons(count).NOTE_clust_qual_confirmed=''; % cluster_quality_confirmed
        NeuronDatabase.neurons(count).NOTE_all_songs_gotten=''; % all_songs_gotten
        NeuronDatabase.neurons(count).NOTE_random=''; % random note
        NeuronDatabase.neurons(count).LEARN_WNonDatestr= tmpstruct.LEARN_WNonDatestr;
        NeuronDatabase.neurons(count).LEARN_WNotherImportantDates=tmpstruct.LEARN_WNotherImportantDates; % leave empty if nothing.
        NeuronDatabase.neurons(count).NOTE_all_labeled=''; %
        
        count = count +1;
    end
end






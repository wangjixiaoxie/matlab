function lt_neural_v2_RemoveNeuron
%% LT, use to manually remove neurons from summary struct

%% provide list to use
SummaryStruct = lt_neural_v2_LoadSummary;



%% remove

birdlist = [];
neuronlist = [];

numbirds = length(SummaryStruct.birds);
count = 1;
for i=1:numbirds
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for j=1:numneurons    
        
        dirname = SummaryStruct.birds(i).neurons(j).dirname;
        disp([num2str(count) ' - ' dirname]);
        
        % ---
        birdlist = [birdlist i];
        neuronlist = [neuronlist j];
        count = count +1;
    end 
end

neuronToRemove = input('HEY! which neuron to remove? (e.g 1)');

%% remove

birdnum = birdlist(neuronToRemove)
neuronnum = neuronlist(neuronToRemove)

SummaryStruct.birds(birdnum).neurons(neuronnum) = [];
disp(['REMOVED neuron ' num2str(neuronToRemove)]);

%% save

targdir = ['/bluejay5/lucas/analyses/neural/'];
save([targdir 'SummaryStruct'], 'SummaryStruct');
disp('DONE!! saved...')


disp('SAVED!');


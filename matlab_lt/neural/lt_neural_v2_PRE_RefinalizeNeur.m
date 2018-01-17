%% currently just a script requiring modificatino for general use.
% refinalize specific neurons

birdnum=1;
neuronlist=[1:8]; % neuron num

for i = neuronlist;
    cd(SummaryStruct.birds(birdnum).neurons(i).dirname);
    %         SummaryStruct.birds(birdnum).neurons(i);
    
    clustnum = SummaryStruct.birds(birdnum).neurons(i).clustnum;
    depth = SummaryStruct.birds(birdnum).neurons(i).electrode_depth;
    Notes = SummaryStruct.birds(birdnum).neurons(i).Notes;
    LEARN_WNonDatestr = SummaryStruct.birds(birdnum).neurons(i).LEARN_WNonDatestr;
    LEARN_WNotherImportantDates = SummaryStruct.birds(birdnum).neurons(i).LEARN_WNotherImportantDates;
    
    
    % CHANGES
    Notes{2} = 'Location_LMAN';
    %         LEARN_WNonDatestr = '15Mar2017-1427';
    
    % RUN
    lt_neural_v2_Finalize(clustnum, depth, Notes, {LEARN_WNonDatestr, LEARN_WNotherImportantDates})
end

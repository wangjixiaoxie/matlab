function ALLBRANCH = lt_neural_v2_CTXT_BranchRemvOlap(ALLBRANCH)

%% lt 9/27/17 - removes overlapping neurons 
% (i.e. if 2 neurons with same data, then removes data for neuron with less
% data

% updates both data and Summarystruct

%% get neurosn to remove and updated summarystruct

[SummaryStruct, NeurToRemove] = lt_neural_v2_DIAGN_RemoveOlap(ALLBRANCH.SummaryStruct);

ALLBRANCH.SummaryStruct = SummaryStruct;


%%  got thru all analysis and remove neurons

numalign = length(ALLBRANCH.alignpos);

for i=1:numalign
    numbirds= length(ALLBRANCH.alignpos(i).bird);
    
    assert(length(NeurToRemove)==numbirds, 'bird missing? bird indices might not match')
    
    for ii=1:numbirds
        
        if isempty(NeurToRemove{ii})
            % then no neurons to remove, skip
            continue
        end
       
       % ====== go thru all branches. for each branch remove appropriate
       % neurons
       numbranches = length(ALLBRANCH.alignpos(i).bird(ii).branch);
       for iii=1:numbranches
           
           if isempty(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron)
               continue
           end
           
           disp(['removed neurons ' NeurToRemove{ii}]);
           
           nremove = NeurToRemove{ii}(NeurToRemove{ii} <= ...
               length(ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron)); % since some neurons might not be present for this branch.
           
           ALLBRANCH.alignpos(i).bird(ii).branch(iii).neuron(nremove) = [];
           
           
       end
    end
end

disp('DONE!!!')
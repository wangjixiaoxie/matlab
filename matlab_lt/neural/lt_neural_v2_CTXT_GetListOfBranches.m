function ALLBRANCH = lt_neural_v2_CTXT_GetListOfBranches(ALLBRANCH)
%% lt 11/29/17 - get lists of regexp strings that make up each branch

% OOUTPUT
% ALLBRANCH.alignpos(apos).bird(i).ListOfBranches

%%
apos = 1;

%%

numbirds = length(ALLBRANCH.alignpos(apos).bird);
for i=1:numbirds
   
    % ============= get list of all branches 
    AllBranches = {};
    numbranches = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for j=1:numbranches 
        if isempty(ALLBRANCH.alignpos(apos).bird(i).branch(j).neuron)
            continue
        end
      
        branches = [ALLBRANCH.alignpos(apos).bird(i).branch(j).neuron.prms_regexpstr];
        AllBranches = [AllBranches branches];
    end

    AllBranches = unique(AllBranches);
    
    ALLBRANCH.alignpos(apos).bird(i).ListOfBranches = AllBranches;
end



function DatStructByBranch = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, ...
    analyfname, tbin)
%%  lt 11/14/17 - if also want to extract premotor window decoding (shuffle), then need:



%% lt 11/29/17 - organizes data by Unique branches (defined by regexpr str)
% OUTPUTS A NEW STRUCTURE


%%
apos = 1; % assume is new analysis version, this always 1.

%% extract all branches first

ALLBRANCH = lt_neural_v2_CTXT_GetListOfBranches(ALLBRANCH);


%% PRODUCE NEW STRUCTURE ORGANIZED BY UNIQUE REGEXP STR

numbirds = length(ALLBRANCH.alignpos(apos).bird);
DatStructByBranch = struct;

for i=1:numbirds
    
    ListOfBranches = ALLBRANCH.alignpos(apos).bird(i).ListOfBranches;
    
    
    % ################################### initiate structure
    for kk=1:length(ListOfBranches)
        DatStructByBranch.bird(i).branchID(kk).regexpstr = '';
        DatStructByBranch.bird(i).branchID(kk).DAT.xdecode = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode_neg = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.ydecode_pos = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.sylcontour_mean = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.sylcontour_x = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.brainregion = {};
        DatStructByBranch.bird(i).branchID(kk).DAT.PREMOTORDECODE_pval = [];
        DatStructByBranch.bird(i).branchID(kk).DAT.PREMOTORDECODE_struct = [];

            DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_apos = [];
            DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_bird = [];
            DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_branch = [];
            DatStructByBranch.bird(i).branchID(kk).DAT.IndsOriginal_neuron = [];

    end
    
    
    % ################################### COLLECT DATA ND PUT INTO STRUCTURE
    numbranches = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for bb =1:numbranches
        
        numneurons = length(ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron);
        for nn=1:numneurons
            
            % ================ collect data
            x_decode = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).xtimes;
            y_decode = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals;
            thisbranch = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).prms_regexpstr;
            
            y_decode_neg = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_neg;
            y_decode_pos = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).yvals_pos;

            sylcontour_mean = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).sylcontours_mean;
            sylcontour_x = ALLBRANCH.alignpos(apos).bird(i).branch(bb).neuron(nn).sylcontours_x;
            
            % ------------------- premotor window shuffles
            if exist('analyfname', 'var')
            decodestruct = lt_neural_v2_CTXT_BRANCH_GetPremotor(analyfname, i, nn, bb, tbin);
            else
            decodestruct = [];    
            end
            if isempty(decodestruct)
                % --- fill with nan
               decodestruct.Pdat = nan;
            end
            
            
            % -------------------- keep data?
            if isempty(x_decode)
                continue
            end
            
            
            
            % ============== save
            ind = strcmp(ListOfBranches, thisbranch);
            DatStructByBranch.bird(i).branchID(ind).regexpstr = thisbranch;
            
            DatStructByBranch.bird(i).branchID(ind).DAT.xdecode = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.xdecode; ...
                x_decode];
           
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode; ...
                y_decode];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_neg= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_neg; ...
                y_decode_neg];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_pos= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.ydecode_pos; ...
                y_decode_pos];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_mean= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_mean; ...
                sylcontour_mean];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_x= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.sylcontour_x; ...
                sylcontour_x];
            
            location = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).NOTE_Location;
            DatStructByBranch.bird(i).branchID(ind).DAT.brainregion = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.brainregion; ...
                location];
            
            % ------------ premotor decode
            DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_pval= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_pval; ...
                decodestruct.Pdat];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_struct= ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.PREMOTORDECODE_struct; ...
                decodestruct];
            
            % ----------- inds to refer back to ALLBRANCH, SUMMARYstruct,
            % ...
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_apos = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_apos; ...
                apos];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_bird = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_bird; ...
                i];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_branch = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_branch; ...
                bb];
            
            DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_neuron = ...
                [DatStructByBranch.bird(i).branchID(ind).DAT.IndsOriginal_neuron; ...
                nn];
            
        end
    end
end
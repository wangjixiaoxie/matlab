function decodestruct = lt_neural_v2_CTXT_BRANCH_GetPremotor(analyfname, birdnum, neurnum, ...
    branchnum, tbin)
%% params
% branchnum; % NOTE: this is not branchID (ie. diff strings can be in same branchnum - 
% this ist he branchnum used in the struct ALLBRANCH

savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

% ========= OUT
% decodestruct = []; if file doesnt exist, otherwise is data.

%% lt 12/14/17 - load previously done shuffle of decode for premotor activity

% ============= filename is:
            fname = [savedir '/' analyfname '/SHUFFDECODE/' ...
                'bird' num2str(birdnum) '_neur' num2str(neurnum) '_branch' num2str(branchnum) ...
                '_tbin' num2str(tbin) '.mat'];
            
            % ============= load
            if ~exist(fname, 'file')
            disp('skipping - no decode dat');
            decodestruct = [];
            else
                disp('extraction succesfful (decode struct)')
            decodestruct = load(fname);
            decodestruct = decodestruct.decodestruct;
            end
            

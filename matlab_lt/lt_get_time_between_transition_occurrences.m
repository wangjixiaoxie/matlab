% LT 9/30/13 - derived from getlabels.m
function [onsets_array offsets_array] =lt_get_time_between_transition_occurrences(bt)

    files=load_batchf(bt);
    
    for ii=1:length(files)
        try
            fn=files(ii).name;
            
            matfn=[fn '.not.mat'];
            strcmd2=['load ' matfn ';']
            eval(strcmd2)
            onsets_array{ii}=onsets;
            offsets_array{ii}=offsets;
        catch err
            continue
        end
    end
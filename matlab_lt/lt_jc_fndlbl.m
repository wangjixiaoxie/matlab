%% LT 10/22/14 - modified from jc_fndlbl to replace each individual instance of a false postiive with a specified replacement syllable.  Each one can be different.

function lt_jc_fndlbl(vlsorind, vlsorfn, ind, newsyll, one_or_multiple_syls)

% LT: added
% one_or_multiple_syls = 0 (use one syl to replace all) or 1 (or use potentially different one for each. Those should then be a cell array with each cell holding the correct replacement syl')
% the cell array can be easily taken directly from "labels" in notmat file

for i = 1:length(ind)
    if ~exist([vlsorfn{ind(i)},'.not.mat'])
        disp([vlsorfn{ind(i)} 'does not exist... error'])
    else
        load([vlsorfn{ind(i)},'.not.mat']);
        if one_or_multiple_syls==0;
            labels(vlsorind(ind(i))) = newsyll;
        elseif one_or_multiple_syls==1;
            labels(vlsorind(ind(i))) = newsyll(ind(i));
        end
        save([vlsorfn{ind(i)},'.not.mat'],'labels','-append');
    end
end

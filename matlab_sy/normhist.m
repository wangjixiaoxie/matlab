function [histout] = normhist(hist)
    for ii=1:length(hist)
        for jj=1:length(hist{ii})
            histout{ii}{jj}=hist{ii}{jj}/sum(hist{ii}{jj})
        end
    end
function jc_fndlbl(vlsorind, vlsorfn, ind, newsyll)


for i = 1:length(ind)
    if ~exist([vlsorfn{ind(i)},'.not.mat'])
        disp(vlsorfn{ind(i)})
    else
    load([vlsorfn{ind(i)},'.not.mat']);
    labels(vlsorind(ind(i))) = newsyll;
    save([vlsorfn{ind(i)},'.not.mat'],'labels','-append');
    end
end

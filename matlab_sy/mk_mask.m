function mask=mk_mask(dat,fs,onset,offset,wid);
%mask=mk_mask(dat,fs,onset,offset,wid);
mask=zeros(size(dat));
t = [1:length(dat)].'/fs;
mask = (0.5*(tanh((t-onset)/wid)+1)).*(1.0-0.5*(tanh((t-offset)/wid)+1));

return;
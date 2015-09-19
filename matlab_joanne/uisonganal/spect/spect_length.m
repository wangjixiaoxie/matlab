function [length, nspect]=spect_length(spect, back)


filt=[1,1,1,1,1];
spect=conv(filt, spect);
back=conv(filt, back);
spect=spect-back;
spect=-1*spect;
nspect=spect/max(spect);
nspect1=nspect(25:125);
nspect2=nspect(26:126);
length=sum(abs(nspect1-nspect2));





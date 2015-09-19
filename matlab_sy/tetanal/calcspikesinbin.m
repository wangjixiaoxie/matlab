function [spksinbin]=tetanal(vec,bintimes,fs);
% [data,fs,spkind,spkamp]=tetanal(fname,TH,song_chan,tet_chans);
% 10/05 changed by TW to take peaktopeak amplitude

ind=find(vec>bintimes(1)*fs&vec<bintimes(2)*fs);
spksinbin=length(ind);
return;

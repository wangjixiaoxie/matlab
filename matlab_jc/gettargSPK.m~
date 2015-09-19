function [toffset]=gettargSPK(batch,note,pretimems)

fs=32000;
if isempty(note)
    note='a';
end

[vals,trigs]=triglabel(batch,note,1,1,0,0);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

toffset=(pretimems*(fs/1000))+((toff/1000)*(fs));

function [toffset]=gettarg(batch,note)

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

toffset=240+((toff/1000)*(fs)-512)/4;

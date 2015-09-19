function [toffset]=gettarg2(batch)

fs=32000;


[vals,trigs]=triglabel(batch,'a',1,1,1,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

toffset=240+((toff/1000)*(fs)-512)/4;

function [fvals,FFstats,toff]=summary_statsFF(templa,batch,batchnote,cntrng)
fs=32000;
fvals=findwnote4(batchnote,'a','','',0,[2000 2700],8500,1,'obs0',0);
mk_tempf(batch,templa,2,'obs0');
get_trigt2(batch,cntrng,0.3,128,1,1);

[vals,trigs]=triglabel(batch,'a',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

vvals=getvals(fvals,1,'TRIG');
FFstats=vvals(:,2);
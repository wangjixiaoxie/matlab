function [FFstats]=fast_stats(templa,batch,batchnote,cntrng)
fvals=findwnote4(batchnote,'a','','',0,[2000 2700],8500,1,'obs0',1);
mk_tempf(batch,templa,2,'obs0');
get_trigt2(batch,cntrng,0.3,128,1,1);

[vals,trigs]=triglabel(batch,'a',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end

vvals=getvals(fvals,1,'TRIG');
ind=find(vvals(:,3)==0);
ind2=find(vvals(:,3)~=0);
FFstats=vvals(ind2,2);


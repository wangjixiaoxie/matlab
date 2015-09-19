function [FFstats,avg]=summary_stats2(templa,fvals,batch,cntrng)
% Do this first: fvals=findwnote4(batchnote,'a','','',0,[2000 2700],8500,1,'obs0',ADDX);
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


count=0;
%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    if isempty(find(ind==i))
        count=count+1;
        shifted(count,:)=fvals(count).datt;
    end
end
[pitch,avg]=jc_pitchmat1024(shifted(1:20,:),1024,1020,1,1950,2600,[1 ],'obs0',1);
g=7;


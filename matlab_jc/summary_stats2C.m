function [fvals,FFstats,avg,toffset,pitch]=summary_stats2C(templa,batch,batchnote,cntrng)
fs=32000;
mk_tempf(batch,templa,2,'obs0');
get_trigt2(batch,cntrng,0.3,128,1,1);

fvals=findwnoteJC(batchnote,'a','','',0,[2000 2700],8500,0,'obs0',1);

[vals,trigs]=triglabel(batch,'a',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

vvals=getvals(fvals,1,'');
FFstats=vvals(:,2);


%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,2,2000,2700,[1 2 3],'obs0',1);
toffset=240+((toff/1000)*(fs)-512)/4;
figure;hold on;
plot(avg);
plot(toffset-50,mean(FFstats),'*')
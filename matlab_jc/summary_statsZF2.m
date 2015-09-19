function [shifted,FFstats,avg,toffset,pitch]=summary_statsZF2(templa,batch,batchnote,cntrng)
fs=32000;
fvals=findwnoteJC2(batchnote,'a','','',0,[3500 4200],8500,1,'obs0',0);
mk_tempf(batch,templa,2,'obs0');
get_trigt2(batch,cntrng,0.3,128,1,1);

[vals,trigs]=triglabel(batch,'b',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

vvals=getvals(fvals,1,'TRIG');
FFstats=vvals(:,2);


%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
[pitch]=jc_pitchmat1024(shifted,1024,1020,2,4400,5100,[1],'obs0',1);
avg=mean(pitch');
toffset=((toff/1000)*(fs)-512)/4;
figure;hold on;
plot(avg);
plot(toffset-50,mean(FFstats),'*')
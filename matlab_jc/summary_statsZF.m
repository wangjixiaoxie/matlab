function [fvals,FFstats,avg,toffset,pitch]=summary_statsZF(templa,batch,batchnote,cntrng)
fs=32000;
fvals=findwnoteJC(batchnote,'a','','',0,[3500 4200],8500,1,'obs0',0);
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


%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
[pitch1]=jc_pitchmat1024(shifted(:,1:3300),1024,1020,2,3300,3900,[1],'obs0',1);
[pitch2]=jc_pitchmat1024(shifted(:,2276:8500),1024,1020,2,3500,4200,[1],'obs0',1);
pitch=[pitch1;pitch2];
avg=mean(pitch');
toffset=((toff/1000)*(fs)-512)/4;
figure;hold on;
plot(avg);
plot(toffset-50,mean(FFstats),'*')
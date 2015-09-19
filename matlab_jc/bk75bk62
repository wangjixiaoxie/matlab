function [fvals,FFstats,avg,toffset,pitch]=summary_statsJC(templa,batch,batchnote,cntrng)
fs=32000;
fvals=findwnoteJC(batchnote,'a','','',0,[2000 2700],8500,1,'obs0',1);
mk_tempf(batch,templa,2,'obs0');
get_trigt2(batch,cntrng,0.2,128,1,1);

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
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,6500,8100,[1],'obs0',1);
% 512 pts - windowing of jc_pitchmat1024
% 240 from 30ms in findwnoteJC
toffset=240+((toff/1000)*(fs)-512)/4;
figure;hold on;
plot(avg);
plot(toffset,mean(FFstats),'*')
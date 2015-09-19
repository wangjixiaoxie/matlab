function [fvals,FFstats,toffset,pitch,avg]=summary_statsTW1(batch,batchnote)
fs=32000;
fvals=findwnoteJC(batchnote,'l','','',0,[3200 3800],8500,1,'obs0',0);

[vals,trigs]=triglabel(batch,'l',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

vvals=getvals(fvals,1,'TRIG');
FFstats=vvals(:,2);

%toffset=((toff/1000)*(fs)-512)/4; %add 160 b/c use findwnoteJC
%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,3200,3800,[3],'obs0',1);
toffset=((toff/1000)*(fs)-512)/4+240;
figure;hold on;
plot(avg);
plot(toffset,mean(FFstats),'*')
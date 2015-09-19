% bk63
fvals720Ahit=findwnoteJC('batch720Anotes','a','','',0,[2000 2700],8500,1,'obs0');
fvals720Amiss=findwnoteJC('batch720Anotes','b','','',0,[2000 2700],8500,1,'obs0');
for i=1:length(fvals720Ahit)
shifted720Ahit(i,:)=fvals720Ahit(i).datt;
end
for i=1:length(fvals720Amiss)
shifted720Amiss(i,:)=fvals720Amiss(i).datt;
end
pitch720Ahit=jc_pitchmat1024(shifted720Ahit,1024,1020,1,2000,2700,[1],'obs0',1);
pitch720Amiss=jc_pitchmat1024(shifted720Amiss,1024,1020,1,2000,2700,[1],'obs0',1);
pitch720A=[pitch720Ahit,pitch720Amiss];
figure;plot(mean(pitchBase713'))
hold on;plot(mean(pitch720A'),'r')
hold on;plot(mean(pitch720Amiss'),'g')
hold on;plot(mean(pitch720Ahit'),'k')

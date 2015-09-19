% bk80
fvals807Ahit=findwnoteJC('batch807Anotes','a','','',0,[2000 2700],8500,1,'obs0');
fvals807Amiss=findwnoteJC('batch807Anotes','b','','',0,[2000 2700],8500,1,'obs0');
for i=1:length(fvals807Ahit)
shifted807Ahit(i,:)=fvals807Ahit(i).datt;
end
for i=1:length(fvals807Amiss)
shifted807Amiss(i,:)=fvals807Amiss(i).datt;
end
pitch807Ahit=jc_pitchmat1024(shifted807Ahit,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807Amiss=jc_pitchmat1024(shifted807Amiss,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807A=[pitch807Ahit,pitch807Amiss];
figure;plot(mean(pitchBase'))
hold on;plot(mean(pitch807A'),'r')
hold on;plot(mean(pitch807Amiss'),'g')
hold on;plot(mean(pitch807Ahit'),'k')

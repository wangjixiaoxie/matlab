mbase1=mean(mean(DShifts(8).pitchBaseline(270:340,:)));
mbase2=mean(mean(DShifts(8).pitchBaseline(462:532,:)));
indEsc=find(mean(DShifts(8).pitchBaseline(270:340,:))<mbase1 & mean(DShifts(8).pitchBaseline(462:532,:))>mbase2);
indHit=find(mean(DShifts(8).pitchBaseline(270:340,:))>mbase1 | mean(DShifts(8).pitchBaseline(462:532,:))<mbase2);

mnbase=mean(DShifts(8).pitchBaseline');
resbase=jc_residuals(DShifts(8).pitchBaseline);
mnlearned=mean(DShifts(8).pitchALL(:,end-200:end)')/mnbase-1;




figure;plot(resbase(:,indHit),'r','Linewidth',2)
hold on;plot(resbase(:,indEsc),'g','Linewidth',2)
mnlearned=mean(DShifts(8).pitchALL(:,end-400:end-200)')./mnbase-1;
hold on;plot(mnlearned,'k','Linewidth',4)
xlim([240 650])
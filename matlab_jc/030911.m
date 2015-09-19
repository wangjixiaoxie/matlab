figure;hold on;
subplot(121);plot([1/8:1/8:161/8],Alldata(1).pitchUDpre(190:350,21:35),'Linewidth',2);ylim([2350 2570]);xlim([0 20])
subplot(122);plot([1/8:1/8:161/8],Alldata(1).pitchUDpost(190:350,21:35),'Linewidth',2);ylim([2350 2570]);xlim([0 20])

mpre=mean(Alldata(1).pitchUDpre(190:350,:));
mpost=mean(Alldata(1).pitchUDpost(190:350,:));
[b,a]=hist(mpre,[2340:20:2580]);
[c,d]=hist(mpost,[2340:20:2580]);
figure;hold on;
stairs(b,a,'Linewidth',2)
stairs(c,d,'r','Linewidth',2)
function jc_plot522(PC,scoresPRE,scoresPOST)
mpr1=mean(scoresPRE(:,1));
mpr2=mean(scoresPRE(:,2));
mpr3=mean(scoresPRE(:,3));
mpo1=mean(scoresPOST(:,1));
mpo2=mean(scoresPOST(:,2));
mpo3=mean(scoresPOST(:,3));
spr1=std(scoresPRE(:,1));
spr2=std(scoresPRE(:,2));
spr3=std(scoresPRE(:,3));
spo1=std(scoresPOST(:,1));
spo2=std(scoresPOST(:,2));
spo3=std(scoresPOST(:,3));
pre=PC(:,1)*mpr1+PC(:,2)*mpr2+PC(:,3)*mpr3;
prePLUS=PC(:,1)*(mpr1+spr1)+PC(:,2)*(mpr2+spr2)+PC(:,3)*(mpr3+spr3);
preMINUS=PC(:,1)*(mpr1-spr1)+PC(:,2)*(mpr2-spr2)+PC(:,3)*(mpr3-spr3);
post=PC(:,1)*mpo1+PC(:,2)*mpo2+PC(:,3)*mpo3;
postPLUS=PC(:,1)*(mpo1+spo1)+PC(:,2)*(mpo2+spo2)+PC(:,3)*(mpo3+spo3);
postMINUS=PC(:,1)*(mpo1-spo1)+PC(:,2)*(mpo2-spo2)+PC(:,3)*(mpo3-spo3);
figure; hold on
plot(pre,'Color','red')
plot(post,'Color','blue')
plot(prePLUS,'Color','red')
plot(preMINUS,'Color','red')
plot(postPLUS,'Color','blue')
plot(postMINUS,'Color','blue')
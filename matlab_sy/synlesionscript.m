%used for making last minute bar plots
%for syntax data
%for sfn
%10.14.09

%synscript
clear mnout
%xvalsprelesion
xprelesion=[125 105 80 80 80 75 75 70 70 60 60 55 30 30]
xpstlesion=[35 20 20 20 20 20 20 20 -10 -10 -10 -10 05 0]


prelesionstartsyn=[.65 .65 .85 .35 .45 .65 .7 .6]
prelesionendsyn=[.05 .25 .05 .2  .25 .05 .05 .35]

pstlesionstartsyn=[.65 .65 .4 .6]
pstlesionendsyn=[.5 .5 .35 .6]

preinds=[1 1 1 8]


wnrawpstlesion=[.15 .15 0 .04] 

%for lesion figure
figure
mnout(1)=mean(xprelesion);
mnout(2)=mean(xpstlesion);
xvls=[1 3]
bar(xvls,mnout,'FaceColor','none','Linewidth',2);
hold on;
%vertical error bar lines
sterr(1)=std(xprelesion)./sqrt(length(xprelesion));
sterr(2)=std(xpstlesion)./sqrt(length(xpstlesion));
for ii=1:2
    plot([xvls(ii) xvls(ii)],[mnout(ii)-sterr(ii) mnout(ii)+sterr(ii)],'Linewidth',2,'Color','k')
end

figure
pctdiffpre=[(prelesionstartsyn-prelesionendsyn)./prelesionstartsyn]
pctdiffpst=[(pstlesionstartsyn-pstlesionendsyn)./pstlesionstartsyn]
pctdiffpre=pctdiffpre*100;
pctdiffpst=pctdiffpst*100
mnoutsyn(1)=mean(pctdiffpre);
mnoutsyn(2)=mean(pctdiffpst);
sterrsyn(1)=std(pctdiffpre)./sqrt(length(pctdiffpre));
sterrsyn(2)=std(pctdiffpst)./sqrt(length(pctdiffpst));
bar(xvls,mnoutsyn,'FaceColor','none','Linewidth',2);
hold on;
for ii=1:2
    plot([xvls(ii) xvls(ii)],[mnoutsyn(ii)-sterrsyn(ii) mnoutsyn(ii)+sterrsyn(ii)],'Linewidth',2,'Color','k')
end



% 
% 
% figure
% pctdiffpre=[xrawprelesion - wnrawprelesion]
% pctdiffpst=[xrawprelesion(preinds)-wnrawpstlesion]
% pctdiffpre=pctdiffpre*100;
% wnrawpstlesion=wnrawpstlesion*100
% mnoutsyn(1)=mean(pctdiffpre);
% mnoutsyn(2)=mean(wnrawpstlesion);
% sterrsyn(1)=std(pctdiffpre)./sqrt(length(pctdiffpre));
% sterrsyn(2)=std(wnrawpstlesion)./sqrt(length(wnrawpstlesion));
% bar(xvls,mnoutsyn,'FaceColor','none','Linewidth',2);
% hold on;
% for ii=1:2
%     plot([xvls(ii) xvls(ii)],[mnoutsyn(ii)-sterrsyn(ii) mnoutsyn(ii)+sterrsyn(ii)],'Linewidth',2,'Color','k')
% end
% 
%         
% % plot([xvals(ii) xvals(ii)], [plotmean(ii)-plotstd(ii) plotmean(ii)+plotstd(ii)],'k');
figure;plot(std(BFLesion(1).UDpre20_1024'))
hold on;plot(std(BFLesion(1).UDpost20_1024'),'r')
for j=1:length(BFLesion)
    for i=1:size(BFLesion(j).UDpre20_1024,1)-280
        centeringBF(j).cvvals(i,:)=mean(std(BFLesion(j).UDpost20_1024(i:i+279,:)'))+mean(std(BFLesion(j).UDpre20_1024(i:i+279,:)'));
    end
    [peak,ind1BF(j)]=min(centeringBF(j).cvvals);
    ind2BF(j)=ind1BF(j)+279;
end

% Determine if the CV threshold is met - 20% reduction - all but 5 and 8
for j=1:length(BFLesion)
    CVpreBF(j)=mean(std(BFLesion(j).UDpre20_1024(ind1BF(j):ind2BF(j),:)'));
    CVpostBF(j)=mean(std(BFLesion(j).UDpost20_1024(ind1BF(j):ind2BF(j),:)'));
end
ratioBF=CVpostBF./CVpreBF;
successesBF=find(ratioBF<0.8);

% Over this window, do it up
for k=1:length(successesBF)
    i=successesBF(k);
    [xB,psdsPREBF(i).data]=jcpsd3(AlldataBFlesion(i).pitchUDpre(ind1BF(i):ind2BF(i),:),8000);
    [xB,psdsPOSTBF(i).data]=jcpsd3(AlldataBFlesion(i).pitchUDpost(ind1BF(i):ind2BF(i),:),8000);
end

for k=1:length(successesBF)
    i=successesBF(k);
    [xB20,psdsPREBF20ms(i).data]=jcpsd3(AlldataBFlesion(i).pitchUDpre(ind1BF(i)+60:ind1BF(i)+219,:),8000);
    [xB20,psdsPOSTBF20ms(i).data]=jcpsd3(AlldataBFlesion(i).pitchUDpost(ind1BF(i)+60:ind1BF(i)+219,:),8000);
end

figure;hold on
for i=1:length(psdsPREBF)
    plot(xB',median(psdsPREBF(i).data)./max(median(psdsPREBF(i).data)),'k')
    plot(xB',median(psdsPOSTBF(i).data)./max(median(psdsPREBF(i).data)),'r')
end
% for i=1:length(psdsPREBF)
%     plot(xB',median(psdsPREBF(i).data),'*')
%     plot(xB',median(psdsPOSTBF(i).data),'*')
% end
figure;hold on
for i=1:length(psdsPREBF)
    plot(xB20',median(psdsPREBF20ms(i).data)./max(median(psdsPREBF20ms(i).data)),'k')
    plot(xB20',median(psdsPOSTBF20ms(i).data)./max(median(psdsPREBF20ms(i).data)),'r')
end
% for i=1:length(psdsPREBF)
%     plot(xB20',median(psdsPREBF20ms(i).data),'*')
%     plot(xB20',median(psdsPOSTBF20ms(i).data),'*')
% end
% 
% figure;hold on
% for i=1:length(psdsPREBF)
% plot(xB',median(psdsPREBF(i).data)./median(psdsPOSTBF(i).data))
% end


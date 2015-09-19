% Determine the 35ms segment in the middle of the note
figure;plot(std(BFLesion(1).UDpre20_1024'))
hold on;plot(std(BFLesion(1).UDpost20_1024'),'r')
for j=1:length(BFLesion)
    for i=1:size(BFLesion(j).UDpre20_1024,1)-160
        centeringBF(j).cvvals(i,:)=mean(std(BFLesion(j).UDpost20_1024(i:i+159,:)'))+mean(std(BFLesion(j).UDpre20_1024(i:i+159,:)'));
    end
    [peak,ind1BF(j)]=min(centeringBF(j).cvvals);
    ind2BF(j)=ind1BF(j)+159;
end

% Determine if the CV threshold is met - 10% reduction - all but 5 and 8
for j=1:length(BFLesion)
    CVpreBF(j)=mean(std(BFLesion(j).UDpre10_512(ind1BF(j):ind2BF(j),:)'));
    CVpostBF(j)=mean(std(BFLesion(j).UDpost10_512(ind1BF(j):ind2BF(j),:)'));
end
ratioBF=CVpostBF./CVpreBF;
successesBF=find(ratioBF<0.8);

% Over this window, do it up
for k=1:length(successesBF)
    i=successesBF(k);
    [xB,psdsPREBF(i).data]=jcpsd2(BFLesion(i).UDpre10_512(ind1BF(i):ind2BF(i),:),8000);
    [xB,psdsPOSTBF(i).data]=jcpsd2(BFLesion(i).UDpost10_512(ind1BF(i):ind2BF(i),:),8000);
end

for k=1:length(successesBF)
    i=successesBF(k);
    [xB10,psdsPREBF10ms(i).data]=jcpsd3(BFLesion(i).UDpre10_512(ind1BF(i)+60:ind1BF(i)+219,:),8000);
    [xB10,psdsPOSTBF10ms(i).data]=jcpsd3(BFLesion(i).UDpost10_512(ind1BF(i)+60:ind1BF(i)+219,:),8000);
end

figure;hold on
for i=1:length(psdsPREBF)
    plot(xB',median(psdsPREBF(i).data),'k')
    plot(xB',median(psdsPOSTBF(i).data),'r')
end
for i=1:length(psdsPREBF)
    plot(xB',median(psdsPREBF(i).data),'*')
    plot(xB',median(psdsPOSTBF(i).data),'*')
end
figure;hold on
for i=1:length(psdsPREBF)
    plot(xB10',median(psdsPREBF10ms(i).data),'k')
    plot(xB10',median(psdsPOSTBF10ms(i).data),'r')
end
for i=1:length(psdsPREBF)
    plot(xB10',median(psdsPREBF10ms(i).data),'*')
    plot(xB10',median(psdsPOSTBF10ms(i).data),'*')
end

figure;hold on
for i=1:length(psdsPRE)
plot(xZ',median(psdsPRE(i).data)./median(psdsPOST(i).data),'r')
end

clear acorr
for i=1:length(successes)
    clear autoc
    for j=1:size(BFLesion(i).UDpost10_512,2)
        autoc(j,:)=xcorr(BFLesion(i).UDpost10_512(:,j));
    end
    acorr(i,:)=median(autoc);
end



















BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,1010,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,508,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,1010,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,508,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,1010,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,508,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,1010,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,508,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');

BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,1010,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,508,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,1010,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,508,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,1010,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,508,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,1010,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,508,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'obs0');

BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,1010,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,508,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,1010,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,508,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,1010,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,508,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,1010,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,508,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'obs0');

BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,1010,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,508,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,1010,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,508,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,1010,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,508,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,1010,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,508,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,1010,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,508,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,1010,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,508,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,1010,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,508,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,1010,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,508,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,1010,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,508,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,1010,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,508,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,1010,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,508,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,1010,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,508,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');

BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,1010,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,508,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,1010,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,508,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,1010,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,508,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,1010,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,508,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,1010,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,508,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,1010,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,508,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,1010,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,508,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,1010,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,508,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');

BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,1010,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,508,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,1010,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,508,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,1010,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,508,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,1010,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,508,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(10).UDpre10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpre,512,1010,1,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpre10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpre,512,508,1,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpre10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpre,512,1010,2,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpre10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpre,512,508,2,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpost10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpost,512,1010,1,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpost10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpost,512,508,1,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpost10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpost,512,1010,2,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');
BFLesion(10).UDpost10_512=jc_PitchData108(AlldataBFlesion(10).rawdataUDpost,512,508,2,AlldataBFlesion(10).FFrange(1),AlldataBFlesion(10).FFrange(2),1,'obs0');

BFLesion(11).UDpre10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpre,512,1010,1,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpre10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpre,512,508,1,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpre10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpre,512,1010,2,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpre10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpre,512,508,2,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpost10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpost,512,1010,1,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpost10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpost,512,508,1,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpost10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpost,512,1010,2,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');
BFLesion(11).UDpost10_512=jc_PitchData108(AlldataBFlesion(11).rawdataUDpost,512,508,2,AlldataBFlesion(11).FFrange(1),AlldataBFlesion(11).FFrange(2),1,'obs0');

BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,1010,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,508,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,1010,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpre,512,508,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,1010,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,508,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,1010,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData108(AlldataBFlesion(1).rawdataUDpost,512,508,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');

BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,1010,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,508,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,1010,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpre10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpre,512,508,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,1010,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,508,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,1010,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost10_512=jc_PitchData108(AlldataBFlesion(2).rawdataUDpost,512,508,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');

BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,1010,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,508,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,1010,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpre10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpre,512,508,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,1010,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,508,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,1010,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost10_512=jc_PitchData108(AlldataBFlesion(3).rawdataUDpost,512,508,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');

BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,1010,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,508,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,1010,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpre10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpre,512,508,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,1010,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,508,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,1010,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost10_512=jc_PitchData108(AlldataBFlesion(4).rawdataUDpost,512,508,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,1010,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,508,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,1010,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpre,512,508,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,1010,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,508,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,1010,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData108(AlldataBFlesion(5).rawdataUDpost,512,508,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,1010,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,508,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,1010,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpre,512,508,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,1010,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,508,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,1010,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData108(AlldataBFlesion(6).rawdataUDpost,512,508,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');

BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,1010,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,508,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,1010,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpre,512,508,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,1010,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,508,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,1010,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData108(AlldataBFlesion(7).rawdataUDpost,512,508,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,1010,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,508,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,1010,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpre,512,508,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,1010,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,508,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,1010,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData108(AlldataBFlesion(8).rawdataUDpost,512,508,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');

BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,1010,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,508,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,1010,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpre,512,508,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,1010,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,508,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,1010,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData108(AlldataBFlesion(9).rawdataUDpost,512,508,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
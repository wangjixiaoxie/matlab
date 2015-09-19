% Determine the 35ms segment in the middle of the note
figure;plot(std(ZFLesion(1).UDpre20_1024'))
hold on;plot(std(ZFLesion(1).UDpost20_1024'),'r')
for j=1:length(ZFLesion)
    for i=1:size(ZFLesion(j).UDpre20_1024,1)-280
        centering(j).cvvals(i,:)=mean(std(ZFLesion(j).UDpost20_1024(i:i+279,:)'))+mean(std(ZFLesion(j).UDpre20_1024(i:i+279,:)'));
    end
    [peak,ind1(j)]=min(centering(j).cvvals);
    ind2(j)=ind1(j)+279;
end

% Determine if the CV threshold is met - 20% reduction - all but 5 and 8
for j=1:length(ZFLesion)
    CVpre(j)=mean(std(ZFLesion(j).UDpre20_1024(ind1(j):ind2(j),:)'));
    CVpost(j)=mean(std(ZFLesion(j).UDpost20_1024(ind1(j):ind2(j),:)'));
end
ratio=CVpost./CVpre;
successes=find(ratio<0.8);

% Over this window, do it up
for k=1:length(successes)
    i=successes(k);
    [xZ,psdsPRE(i).data]=jcpsd2(AlldataZFlesion(i).pitchUDpre(ind1(i):ind2(i),:),8000);
    [xZ,psdsPOST(i).data]=jcpsd2(AlldataZFlesion(i).pitchUDpost(ind1(i):ind2(i),:),8000);
end

for k=1:length(successes)
    i=successes(k);
    [xZ20,psdsPRE20ms(i).data]=jcpsd2(AlldataZFlesion(i).UDpre20_1024(ind1(i)+60:ind1(i)+219,:),8000);
    [xZ20,psdsPOST20ms(i).data]=jcpsd2(AlldataZFlesion(i).UDpost20_1024(ind1(i)+60:ind1(i)+219,:),8000);
end

figure;hold on
for i=1:length(psdsPRE)
    plot(xZ',median(psdsPRE(i).data)./max(median(psdsPRE(i).data)),'k')
    plot(xZ',median(psdsPOST(i).data)./max(median(psdsPRE(i).data)),'r')
end
for i=1:length(psdsPRE)
    plot(xZ',median(psdsPRE(i).data),'*')
    plot(xZ',median(psdsPOST(i).data),'*')
end
figure;hold on
for i=1:length(psdsPRE)
    plot(xZ20',median(psdsPRE20ms(i).data)./max(median(psdsPRE20ms(i).data)),'k')
    plot(xZ20',median(psdsPOST20ms(i).data)./max(median(psdsPRE20ms(i).data)),'r')
end
for i=1:length(psdsPRE)
    plot(xZ20',median(psdsPRE20ms(i).data),'*')
    plot(xZ20',median(psdsPOST20ms(i).data),'*')
end

figure;hold on
for i=1:length(psdsPRE)
plot(xZ',median(psdsPRE(i).data)./median(psdsPOST(i).data),'r')
end

clear acorr
for i=1:length(successes)
    clear autoc

        autoc=xcov(ZFLesion(i).UDpost20_1024(ind1(successes(i)):ind2(successes(i)),:));
    acov(i,:)=median(autoc');
end



















ZFLesion(1).UDpre10_1024=jc_PitchData208(AlldataZFlesion(1).rawdataUDpre,1024,1020,1,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpre10_512=jc_PitchData208(AlldataZFlesion(1).rawdataUDpre,512,508,1,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpre20_1024=jc_PitchData208(AlldataZFlesion(1).rawdataUDpre,1024,1020,2,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpre20_512=jc_PitchData208(AlldataZFlesion(1).rawdataUDpre,512,508,2,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpost10_1024=jc_PitchData208(AlldataZFlesion(1).rawdataUDpost,1024,1020,1,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpost10_512=jc_PitchData208(AlldataZFlesion(1).rawdataUDpost,512,508,1,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpost20_1024=jc_PitchData208(AlldataZFlesion(1).rawdataUDpost,1024,1020,2,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');
ZFLesion(1).UDpost20_512=jc_PitchData208(AlldataZFlesion(1).rawdataUDpost,512,508,2,AlldataZFlesion(1).FFrange(1),AlldataZFlesion(1).FFrange(2),1,'obs0');

ZFLesion(2).UDpre10_1024=jc_PitchData208(AlldataZFlesion(2).rawdataUDpre,1024,1020,1,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpre10_512=jc_PitchData208(AlldataZFlesion(2).rawdataUDpre,512,508,1,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpre20_1024=jc_PitchData208(AlldataZFlesion(2).rawdataUDpre,1024,1020,2,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpre20_512=jc_PitchData208(AlldataZFlesion(2).rawdataUDpre,512,508,2,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpost10_1024=jc_PitchData208(AlldataZFlesion(2).rawdataUDpost,1024,1020,1,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpost10_512=jc_PitchData208(AlldataZFlesion(2).rawdataUDpost,512,508,1,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpost20_1024=jc_PitchData208(AlldataZFlesion(2).rawdataUDpost,1024,1020,2,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');
ZFLesion(2).UDpost20_512=jc_PitchData208(AlldataZFlesion(2).rawdataUDpost,512,508,2,AlldataZFlesion(2).FFrange(1),AlldataZFlesion(2).FFrange(2),1,'obs0');

ZFLesion(3).UDpre10_1024=jc_PitchData208(AlldataZFlesion(3).rawdataUDpre,1024,1020,1,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpre10_512=jc_PitchData208(AlldataZFlesion(3).rawdataUDpre,512,508,1,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpre20_1024=jc_PitchData208(AlldataZFlesion(3).rawdataUDpre,1024,1020,2,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpre20_512=jc_PitchData208(AlldataZFlesion(3).rawdataUDpre,512,508,2,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpost10_1024=jc_PitchData208(AlldataZFlesion(3).rawdataUDpost,1024,1020,1,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpost10_512=jc_PitchData208(AlldataZFlesion(3).rawdataUDpost,512,508,1,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpost20_1024=jc_PitchData208(AlldataZFlesion(3).rawdataUDpost,1024,1020,2,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');
ZFLesion(3).UDpost20_512=jc_PitchData208(AlldataZFlesion(3).rawdataUDpost,512,508,2,AlldataZFlesion(3).FFrange(1),AlldataZFlesion(3).FFrange(2),1,'obs0');

ZFLesion(4).UDpre10_1024=jc_PitchData208(AlldataZFlesion(4).rawdataUDpre,1024,1020,1,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpre10_512=jc_PitchData208(AlldataZFlesion(4).rawdataUDpre,512,508,1,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpre20_1024=jc_PitchData208(AlldataZFlesion(4).rawdataUDpre,1024,1020,2,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpre20_512=jc_PitchData208(AlldataZFlesion(4).rawdataUDpre,512,508,2,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpost10_1024=jc_PitchData208(AlldataZFlesion(4).rawdataUDpost,1024,1020,1,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpost10_512=jc_PitchData208(AlldataZFlesion(4).rawdataUDpost,512,508,1,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpost20_1024=jc_PitchData208(AlldataZFlesion(4).rawdataUDpost,1024,1020,2,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(4).UDpost20_512=jc_PitchData208(AlldataZFlesion(4).rawdataUDpost,512,508,2,AlldataZFlesion(4).FFrange(1),AlldataZFlesion(4).FFrange(2),1,'obs0');
ZFLesion(5).UDpre10_1024=jc_PitchData208(AlldataZFlesion(5).rawdataUDpre,1024,1020,1,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpre10_512=jc_PitchData208(AlldataZFlesion(5).rawdataUDpre,512,508,1,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpre20_1024=jc_PitchData208(AlldataZFlesion(5).rawdataUDpre,1024,1020,2,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpre20_512=jc_PitchData208(AlldataZFlesion(5).rawdataUDpre,512,508,2,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpost10_1024=jc_PitchData208(AlldataZFlesion(5).rawdataUDpost,1024,1020,1,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpost10_512=jc_PitchData208(AlldataZFlesion(5).rawdataUDpost,512,508,1,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpost20_1024=jc_PitchData208(AlldataZFlesion(5).rawdataUDpost,1024,1020,2,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(5).UDpost20_512=jc_PitchData208(AlldataZFlesion(5).rawdataUDpost,512,508,2,AlldataZFlesion(5).FFrange(1),AlldataZFlesion(5).FFrange(2),1,'obs0');
ZFLesion(6).UDpre10_1024=jc_PitchData208(AlldataZFlesion(6).rawdataUDpre,1024,1020,1,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpre10_512=jc_PitchData208(AlldataZFlesion(6).rawdataUDpre,512,508,1,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpre20_1024=jc_PitchData208(AlldataZFlesion(6).rawdataUDpre,1024,1020,2,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpre20_512=jc_PitchData208(AlldataZFlesion(6).rawdataUDpre,512,508,2,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpost10_1024=jc_PitchData208(AlldataZFlesion(6).rawdataUDpost,1024,1020,1,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpost10_512=jc_PitchData208(AlldataZFlesion(6).rawdataUDpost,512,508,1,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpost20_1024=jc_PitchData208(AlldataZFlesion(6).rawdataUDpost,1024,1020,2,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');
ZFLesion(6).UDpost20_512=jc_PitchData208(AlldataZFlesion(6).rawdataUDpost,512,508,2,AlldataZFlesion(6).FFrange(1),AlldataZFlesion(6).FFrange(2),1,'obs0');

ZFLesion(7).UDpre10_1024=jc_PitchData208(AlldataZFlesion(7).rawdataUDpre,1024,1020,1,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpre10_512=jc_PitchData208(AlldataZFlesion(7).rawdataUDpre,512,508,1,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpre20_1024=jc_PitchData208(AlldataZFlesion(7).rawdataUDpre,1024,1020,2,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpre20_512=jc_PitchData208(AlldataZFlesion(7).rawdataUDpre,512,508,2,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpost10_1024=jc_PitchData208(AlldataZFlesion(7).rawdataUDpost,1024,1020,1,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpost10_512=jc_PitchData208(AlldataZFlesion(7).rawdataUDpost,512,508,1,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpost20_1024=jc_PitchData208(AlldataZFlesion(7).rawdataUDpost,1024,1020,2,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(7).UDpost20_512=jc_PitchData208(AlldataZFlesion(7).rawdataUDpost,512,508,2,AlldataZFlesion(7).FFrange(1),AlldataZFlesion(7).FFrange(2),1,'obs0');
ZFLesion(8).UDpre10_1024=jc_PitchData208(AlldataZFlesion(8).rawdataUDpre,1024,1020,1,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpre10_512=jc_PitchData208(AlldataZFlesion(8).rawdataUDpre,512,508,1,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpre20_1024=jc_PitchData208(AlldataZFlesion(8).rawdataUDpre,1024,1020,2,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpre20_512=jc_PitchData208(AlldataZFlesion(8).rawdataUDpre,512,508,2,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpost10_1024=jc_PitchData208(AlldataZFlesion(8).rawdataUDpost,1024,1020,1,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpost10_512=jc_PitchData208(AlldataZFlesion(8).rawdataUDpost,512,508,1,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpost20_1024=jc_PitchData208(AlldataZFlesion(8).rawdataUDpost,1024,1020,2,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');
ZFLesion(8).UDpost20_512=jc_PitchData208(AlldataZFlesion(8).rawdataUDpost,512,508,2,AlldataZFlesion(8).FFrange(1),AlldataZFlesion(8).FFrange(2),1,'obs0');

ZFLesion(9).UDpre10_1024=jc_PitchData208(AlldataZFlesion(9).rawdataUDpre,1024,1020,1,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpre10_512=jc_PitchData208(AlldataZFlesion(9).rawdataUDpre,512,508,1,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpre20_1024=jc_PitchData208(AlldataZFlesion(9).rawdataUDpre,1024,1020,2,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpre20_512=jc_PitchData208(AlldataZFlesion(9).rawdataUDpre,512,508,2,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpost10_1024=jc_PitchData208(AlldataZFlesion(9).rawdataUDpost,1024,1020,1,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpost10_512=jc_PitchData208(AlldataZFlesion(9).rawdataUDpost,512,508,1,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpost20_1024=jc_PitchData208(AlldataZFlesion(9).rawdataUDpost,1024,1020,2,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(9).UDpost20_512=jc_PitchData208(AlldataZFlesion(9).rawdataUDpost,512,508,2,AlldataZFlesion(9).FFrange(1),AlldataZFlesion(9).FFrange(2),1,'obs0');
ZFLesion(10).UDpre10_1024=jc_PitchData208(AlldataZFlesion(10).rawdataUDpre,1024,1020,1,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpre10_512=jc_PitchData208(AlldataZFlesion(10).rawdataUDpre,512,508,1,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpre20_1024=jc_PitchData208(AlldataZFlesion(10).rawdataUDpre,1024,1020,2,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpre20_512=jc_PitchData208(AlldataZFlesion(10).rawdataUDpre,512,508,2,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpost10_1024=jc_PitchData208(AlldataZFlesion(10).rawdataUDpost,1024,1020,1,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpost10_512=jc_PitchData208(AlldataZFlesion(10).rawdataUDpost,512,508,1,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpost20_1024=jc_PitchData208(AlldataZFlesion(10).rawdataUDpost,1024,1020,2,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');
ZFLesion(10).UDpost20_512=jc_PitchData208(AlldataZFlesion(10).rawdataUDpost,512,508,2,AlldataZFlesion(10).FFrange(1),AlldataZFlesion(10).FFrange(2),1,'obs0');

ZFLesion(11).UDpre10_1024=jc_PitchData208(AlldataZFlesion(11).rawdataUDpre,1024,1020,1,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpre10_512=jc_PitchData208(AlldataZFlesion(11).rawdataUDpre,512,508,1,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpre20_1024=jc_PitchData208(AlldataZFlesion(11).rawdataUDpre,1024,1020,2,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpre20_512=jc_PitchData208(AlldataZFlesion(11).rawdataUDpre,512,508,2,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpost10_1024=jc_PitchData208(AlldataZFlesion(11).rawdataUDpost,1024,1020,1,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpost10_512=jc_PitchData208(AlldataZFlesion(11).rawdataUDpost,512,508,1,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpost20_1024=jc_PitchData208(AlldataZFlesion(11).rawdataUDpost,1024,1020,2,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');
ZFLesion(11).UDpost20_512=jc_PitchData208(AlldataZFlesion(11).rawdataUDpost,512,508,2,AlldataZFlesion(11).FFrange(1),AlldataZFlesion(11).FFrange(2),1,'obs0');

BFLesion(1).UDpre10_1024=jc_PitchData208(AlldataBFlesion(1).rawdataUDpre,1024,1020,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre10_512=jc_PitchData208(AlldataBFlesion(1).rawdataUDpre,512,508,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre20_1024=jc_PitchData208(AlldataBFlesion(1).rawdataUDpre,1024,1020,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpre20_512=jc_PitchData208(AlldataBFlesion(1).rawdataUDpre,512,508,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_1024=jc_PitchData208(AlldataBFlesion(1).rawdataUDpost,1024,1020,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost10_512=jc_PitchData208(AlldataBFlesion(1).rawdataUDpost,512,508,1,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost20_1024=jc_PitchData208(AlldataBFlesion(1).rawdataUDpost,1024,1020,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');
BFLesion(1).UDpost20_512=jc_PitchData208(AlldataBFlesion(1).rawdataUDpost,512,508,2,AlldataBFlesion(1).FFrange(1),AlldataBFlesion(1).FFrange(2),1,'obs0');

BFLesion(2).UDpre10_1024=jc_PitchData208(AlldataBFlesion(2).rawdataUDpre,1024,1020,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpre10_512=jc_PitchData208(AlldataBFlesion(2).rawdataUDpre,512,508,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpre20_1024=jc_PitchData208(AlldataBFlesion(2).rawdataUDpre,1024,1020,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpre20_512=jc_PitchData208(AlldataBFlesion(2).rawdataUDpre,512,508,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost10_1024=jc_PitchData208(AlldataBFlesion(2).rawdataUDpost,1024,1020,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost10_512=jc_PitchData208(AlldataBFlesion(2).rawdataUDpost,512,508,1,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost20_1024=jc_PitchData208(AlldataBFlesion(2).rawdataUDpost,1024,1020,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');
BFLesion(2).UDpost20_512=jc_PitchData208(AlldataBFlesion(2).rawdataUDpost,512,508,2,AlldataBFlesion(2).FFrange(1),AlldataBFlesion(2).FFrange(2),1,'w');

BFLesion(3).UDpre10_1024=jc_PitchData208(AlldataBFlesion(3).rawdataUDpre,1024,1020,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpre10_512=jc_PitchData208(AlldataBFlesion(3).rawdataUDpre,512,508,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpre20_1024=jc_PitchData208(AlldataBFlesion(3).rawdataUDpre,1024,1020,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpre20_512=jc_PitchData208(AlldataBFlesion(3).rawdataUDpre,512,508,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost10_1024=jc_PitchData208(AlldataBFlesion(3).rawdataUDpost,1024,1020,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost10_512=jc_PitchData208(AlldataBFlesion(3).rawdataUDpost,512,508,1,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost20_1024=jc_PitchData208(AlldataBFlesion(3).rawdataUDpost,1024,1020,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');
BFLesion(3).UDpost20_512=jc_PitchData208(AlldataBFlesion(3).rawdataUDpost,512,508,2,AlldataBFlesion(3).FFrange(1),AlldataBFlesion(3).FFrange(2),1,'w');

BFLesion(4).UDpre10_1024=jc_PitchData208(AlldataBFlesion(4).rawdataUDpre,1024,1020,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpre10_512=jc_PitchData208(AlldataBFlesion(4).rawdataUDpre,512,508,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpre20_1024=jc_PitchData208(AlldataBFlesion(4).rawdataUDpre,1024,1020,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpre20_512=jc_PitchData208(AlldataBFlesion(4).rawdataUDpre,512,508,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost10_1024=jc_PitchData208(AlldataBFlesion(4).rawdataUDpost,1024,1020,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost10_512=jc_PitchData208(AlldataBFlesion(4).rawdataUDpost,512,508,1,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost20_1024=jc_PitchData208(AlldataBFlesion(4).rawdataUDpost,1024,1020,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(4).UDpost20_512=jc_PitchData208(AlldataBFlesion(4).rawdataUDpost,512,508,2,AlldataBFlesion(4).FFrange(1),AlldataBFlesion(4).FFrange(2),1,'w');
BFLesion(5).UDpre10_1024=jc_PitchData208(AlldataBFlesion(5).rawdataUDpre,1024,1020,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre10_512=jc_PitchData208(AlldataBFlesion(5).rawdataUDpre,512,508,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre20_1024=jc_PitchData208(AlldataBFlesion(5).rawdataUDpre,1024,1020,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpre20_512=jc_PitchData208(AlldataBFlesion(5).rawdataUDpre,512,508,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_1024=jc_PitchData208(AlldataBFlesion(5).rawdataUDpost,1024,1020,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost10_512=jc_PitchData208(AlldataBFlesion(5).rawdataUDpost,512,508,1,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost20_1024=jc_PitchData208(AlldataBFlesion(5).rawdataUDpost,1024,1020,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(5).UDpost20_512=jc_PitchData208(AlldataBFlesion(5).rawdataUDpost,512,508,2,AlldataBFlesion(5).FFrange(1),AlldataBFlesion(5).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_1024=jc_PitchData208(AlldataBFlesion(6).rawdataUDpre,1024,1020,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre10_512=jc_PitchData208(AlldataBFlesion(6).rawdataUDpre,512,508,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre20_1024=jc_PitchData208(AlldataBFlesion(6).rawdataUDpre,1024,1020,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpre20_512=jc_PitchData208(AlldataBFlesion(6).rawdataUDpre,512,508,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_1024=jc_PitchData208(AlldataBFlesion(6).rawdataUDpost,1024,1020,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost10_512=jc_PitchData208(AlldataBFlesion(6).rawdataUDpost,512,508,1,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost20_1024=jc_PitchData208(AlldataBFlesion(6).rawdataUDpost,1024,1020,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');
BFLesion(6).UDpost20_512=jc_PitchData208(AlldataBFlesion(6).rawdataUDpost,512,508,2,AlldataBFlesion(6).FFrange(1),AlldataBFlesion(6).FFrange(2),1,'obs0');

BFLesion(7).UDpre10_1024=jc_PitchData208(AlldataBFlesion(7).rawdataUDpre,1024,1020,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre10_512=jc_PitchData208(AlldataBFlesion(7).rawdataUDpre,512,508,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre20_1024=jc_PitchData208(AlldataBFlesion(7).rawdataUDpre,1024,1020,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpre20_512=jc_PitchData208(AlldataBFlesion(7).rawdataUDpre,512,508,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_1024=jc_PitchData208(AlldataBFlesion(7).rawdataUDpost,1024,1020,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost10_512=jc_PitchData208(AlldataBFlesion(7).rawdataUDpost,512,508,1,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost20_1024=jc_PitchData208(AlldataBFlesion(7).rawdataUDpost,1024,1020,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(7).UDpost20_512=jc_PitchData208(AlldataBFlesion(7).rawdataUDpost,512,508,2,AlldataBFlesion(7).FFrange(1),AlldataBFlesion(7).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_1024=jc_PitchData208(AlldataBFlesion(8).rawdataUDpre,1024,1020,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre10_512=jc_PitchData208(AlldataBFlesion(8).rawdataUDpre,512,508,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre20_1024=jc_PitchData208(AlldataBFlesion(8).rawdataUDpre,1024,1020,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpre20_512=jc_PitchData208(AlldataBFlesion(8).rawdataUDpre,512,508,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_1024=jc_PitchData208(AlldataBFlesion(8).rawdataUDpost,1024,1020,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost10_512=jc_PitchData208(AlldataBFlesion(8).rawdataUDpost,512,508,1,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost20_1024=jc_PitchData208(AlldataBFlesion(8).rawdataUDpost,1024,1020,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');
BFLesion(8).UDpost20_512=jc_PitchData208(AlldataBFlesion(8).rawdataUDpost,512,508,2,AlldataBFlesion(8).FFrange(1),AlldataBFlesion(8).FFrange(2),1,'obs0');

BFLesion(9).UDpre10_1024=jc_PitchData208(AlldataBFlesion(9).rawdataUDpre,1024,1020,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre10_512=jc_PitchData208(AlldataBFlesion(9).rawdataUDpre,512,508,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre20_1024=jc_PitchData208(AlldataBFlesion(9).rawdataUDpre,1024,1020,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpre20_512=jc_PitchData208(AlldataBFlesion(9).rawdataUDpre,512,508,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_1024=jc_PitchData208(AlldataBFlesion(9).rawdataUDpost,1024,1020,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost10_512=jc_PitchData208(AlldataBFlesion(9).rawdataUDpost,512,508,1,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost20_1024=jc_PitchData208(AlldataBFlesion(9).rawdataUDpost,1024,1020,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
BFLesion(9).UDpost20_512=jc_PitchData208(AlldataBFlesion(9).rawdataUDpost,512,508,2,AlldataBFlesion(9).FFrange(1),AlldataBFlesion(9).FFrange(2),1,'obs0');
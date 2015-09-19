
%plot(X,Y,L,U)
ebars1=[std(HistAPV1(end-99:end))/10 std(Histpost1(1:50))/sqrt(50) std(Histpost1(51:100))/sqrt(50)];
ebars2=[std(HistAPV2(end-99:end))/10 std(Histpost2(1:50))/sqrt(50) std(Histpost2(51:100))/sqrt(50)];
figure;errorbar([1.1 2 3],[0 mean(Histpost1(1:50))-mean(HistAPV1(end-99:end)) mean(Histpost1(51:100))-mean(HistAPV1(end-99:end))],ebars,ebars)
hold on;errorbar([0.9 2 3],[0 mean(Histpost2(1:50))-mean(HistAPV2(end-99:end)) mean(Histpost2(51:100))-mean(HistAPV2(end-99:end))],ebars2,ebars2)
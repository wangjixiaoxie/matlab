function plot_errorbars(FFstats,prct)
figure;hold on;

for i=1:length(FFstats)
    GG=FFstats(i).data;
    X(i)=i;
    Y(i)=median(GG);
    L(i)=median(GG)-prctile(GG,prct);
    U(i)=prctile(GG,100-prct)-median(GG);
end

figure;errorbar(X,Y,L,U)
%plot(X,Y,L,U)
load /bulbul3/AnalysisRepeats/Exp_100611.mat
load /bulbul3/AnalysisRepeats/LMANplots
figure;hold on;
subplotindices=[3 8 11 13 16 17 18];
expindices=[1 2 12 14 9 10 15];

for i=1:length(expindices)
    subplot(4,5,subplotindices(i));hold on;
    k=expindices(i);
    plot([Exp1(k).predays Exp1(k).wndays Exp1(k).postdays]-min(Exp1(k).wndays),Exp1(k).replong([Exp1(k).predays Exp1(k).wndays Exp1(k).postdays]),'k')
    plot(Exp1(k).predays-min(Exp1(k).wndays),Exp1(k).replong(Exp1(k).predays),'.','Markersize',15)
    plot(Exp1(k).wndays-min(Exp1(k).wndays),Exp1(k).replong(Exp1(k).wndays),'r.','Markersize',15)
    plot(Exp1(k).postdays-min(Exp1(k).wndays),Exp1(k).replong(Exp1(k).postdays),'k.','Markersize',15)
    xlim(xlimits(i,:))
    ylim(ylimits(i,:))
end

for i=4     % TW has all other data
    subplot(4,5,(i-1)*5+4);hold on;
    plot(timing3(FF(i).prefv{1}),mean(FF(i).prept{1}(FF(i).window,:)),'.')
    plot(timing3(FF(i).prefv{2}),mean(FF(i).prept{2}(FF(i).window,:)),'r.')
    plot(timing3(FF(i).prefv{3}),mean(FF(i).prept{3}(FF(i).window,:)),'k.')
    ylim(FF(i).ylim)    
    subplot(4,5,i*5);hold on;
    plot(timing3(FF(i).postfv{1}),mean(FF(i).postpt{1}(FF(i).window,:)),'.')
    plot(timing3(FF(i).postfv{2}),mean(FF(i).postpt{2}(FF(i).window,:)),'r.')
    plot(timing3(FF(i).postfv{3}),mean(FF(i).postpt{3}(FF(i).window,:)),'k.')
    ylim(FF(i).ylim)
end

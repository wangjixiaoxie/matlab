% 
figure;hold on;
for k=1:length(ind)
    i=ind(k);
    subplot(5,4,k);hold on;
    coef=1-2*(isequal(Experiment(i).DIR,'down'));
    plot(runningaverage(Experiment(i).timeACpre,10),coef*runningaverage(mean(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,:)),10))
    plot(runningaverage(Experiment(i).timeAPV,10),coef*runningaverage(mean(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,:)),10),'g')
    plot(runningaverage(Experiment(i).timeAPVwn,10),coef*runningaverage(mean(Experiment(i).pitchAPVwn(Experiment(i).on:Experiment(i).off,:)),10),'r')
    plot(runningaverage(Experiment(i).timeACpost,10),coef*runningaverage(mean(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,:)),10),'k')
    xlim([mean(Experiment(i).timeAPVwn)-24 mean(Experiment(i).timeAPVwn)+24])
    ylim([coef*mean(mean(Experiment(i).pitchAPVwn(Experiment(i).on:Experiment(i).off,:)))-100 coef*mean(mean(Experiment(i).pitchAPVwn(Experiment(i).on:Experiment(i).off,:)))+100])
end
%  
figure;hold on;
for k=1:length(ind)
    i=ind(k);
    subplot(5,4,k);hold on;
    coef=1-2*(isequal(Experiment(i).DIR,'down'));
    if ~isempty(Experiment(i).timeAPVwnCTL)
    plot(runningaverage(Experiment(i).timeACpreCTL,10),coef*runningaverage(mean(Experiment(i).pitchACpreCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)),10))
    plot(runningaverage(Experiment(i).timeAPVCTL,10),coef*runningaverage(mean(Experiment(i).pitchAPVCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)),10),'g')
    plot(runningaverage(Experiment(i).timeAPVwnCTL,10),coef*runningaverage(mean(Experiment(i).pitchAPVwnCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)),10),'r')
    plot(runningaverage(Experiment(i).timeACpostCTL,10),coef*runningaverage(mean(Experiment(i).pitchACpostCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)),10),'k')
    xlim([mean(Experiment(i).timeAPVwnCTL)-24 mean(Experiment(i).timeAPVwnCTL)+24])
    ylim([coef*mean(mean(Experiment(i).pitchAPVwnCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)))-100 coef*mean(mean(Experiment(i).pitchAPVwnCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)))+100])
    end
    if ~isempty(Experiment(i).timeACpreCTL)
        plot(runningaverage(Experiment(i).timeACpreCTL,10),coef*runningaverage(mean(Experiment(i).pitchACpreCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)),10))
        plot(runningaverage(Experiment(i).timeACpostCTL,10),coef*runningaverage(mean(Experiment(i).pitchACpostCTL(Experiment(i).onCTL:Experiment(i).offCTL,:)),10),'k')
    end
    
    
end
    

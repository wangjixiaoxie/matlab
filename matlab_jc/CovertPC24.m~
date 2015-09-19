% All data
                figure;hold on;
                for j=1:4
                    targmed=round(median(PC24(j).targeting)); % actually 0-64 pts earlier, but syllable is noisy there
                    sigma=std(PC24(j).pitchPRE(targmed,:));
                    s(j)=sigma;
                    subplot(2,2,j);hold on;
                    PREdays=[1 find(diff(PC24(j).tvPRE)>8) length(PC24(j).tvPRE)];
                    WNdays=[1 find(diff(PC24(j).tvWN)>8) length(PC24(j).tvWN)];
                    POSTdays=[1 find(diff(PC24(j).tvPOST)>8) length(PC24(j).tvPOST)];
                    for i=1:length(PREdays)-1
                        if PREdays(i+1)>PREdays(i)+ws
                            plot(runningaverage(PC24(j).tvPRE(PREdays(i):PREdays(i+1))-PC24(j).tvWN(1),ws),coef*runningaverage((PC24(j).pitchPRE(targmed,PREdays(i):PREdays(i+1))-mean(PC24(j).pitchPRE(targmed,:)))/sigma,ws),'*');
                        end
                    end
                    for i=1:length(WNdays)-1
                        if WNdays(i+1)>WNdays(i)+ws
                            plot(runningaverage(PC24(j).tvWN(WNdays(i):WNdays(i+1))-PC24(j).tvWN(1),ws),coef*runningaverage((PC24(j).pitchWN(targmed,WNdays(i):WNdays(i+1))-mean(PC24(j).pitchPRE(targmed,:)))/sigma,ws),'*','Color','r')
                        end
                    end    
                    for i=1:length(POSTdays)-1
                        if POSTdays(i+1)>POSTdays(i)+ws
                            plot(runningaverage(PC24(j).tvPOST(POSTdays(i):POSTdays(i+1))-PC24(j).tvWN(1),ws),coef*runningaverage((PC24(j).pitchPOST(targmed,POSTdays(i):POSTdays(i+1))-mean(PC24(j).pitchPRE(targmed,:)))/sigma,ws),'*','Color','k')
                        end
                    end
                    plot([-25 70],[0 0])
                    thr=0.5*(1-2*isequal(PC24(j).direction,'down'));
                    plot([-25:5:70],[ones(1,20)*thr],'.')
                    xlim([-25 70])
                    ylim([-2 2])
                end


% Targeting, error bars, etc
                for j=1:4
                    targmed=round(median(PC24(j).targeting)); % actually 0-64 pts earlier, but syllable is occasionally noisy there
                    sigma=std(PC24(j).pitchPRE(targmed,:));
                    coef=(1-2*isequal(PC24(j).direction,'down'));
                    firstday=find(PC24(j).tvPOST<(PC24(j).tvPOST(1)+10));
                    avgadaptive(j,:)=coef*(mean(PC24(j).pitchPOST(targmed-32-200:targmed-32+200,firstday)')-mean(PC24(j).pitchPRE(targmed-32-200:targmed-32+200,:)'));
                    thr=0.5*coef;
                end
                figure;hold on;
                for i=1:4
                plot([-25:1/8:25],avgadaptive(i,:))
                end
                %plot([0 0],[0 150],'Color','k','Linewidth',2)
                xlim([-15 25])
                ylim([-20 120])
                errorbar(0,mean(avgadaptive(:,201)),std(avgadaptive(:,201))/2)


% figure;hold on;
% ws=20;
% for j=1:4
%     targmed=round(median(PC24(j).targeting)); % actually 0-64 pts earlier, but syllable is noisy there
%     sigma=std(PC24(j).pitchPRE(targmed,:));
%     coef=(1-2*isequal(PC24(j).direction,'down'));
%     PREdays=[1 find(diff(PC24(j).tvPRE)>8) length(PC24(j).tvPRE)];
%     WNdays=[1 find(diff(PC24(j).tvWN)>8) length(PC24(j).tvWN)];
%     POSTdays=[1 find(diff(PC24(j).tvPOST)>8) length(PC24(j).tvPOST)];
%     for i=1:length(PREdays)-1
%         if PREdays(i+1)>PREdays(i)+ws
%             plot(runningaverage(PC24(j).tvPRE(PREdays(i):PREdays(i+1))-PC24(j).tvWN(1),ws),coef*runningaverage((PC24(j).pitchPRE(targmed,PREdays(i):PREdays(i+1))-mean(PC24(j).pitchPRE(targmed,:)))/sigma,ws));
%         end
%     end
%     for i=1:length(WNdays)-1
%         if WNdays(i+1)>WNdays(i)+ws
%             plot(runningaverage(PC24(j).tvWN(WNdays(i):WNdays(i+1))-PC24(j).tvWN(1),ws),coef*runningaverage((PC24(j).pitchWN(targmed,WNdays(i):WNdays(i+1))-mean(PC24(j).pitchPRE(targmed,:)))/sigma,ws),'Color','r')
%         end
%     end    
%     for i=1:length(POSTdays)-1
%         if POSTdays(i+1)>POSTdays(i)+ws
%             plot(runningaverage(PC24(j).tvPOST(POSTdays(i):POSTdays(i+1))-PC24(j).tvWN(1),ws),coef*runningaverage((PC24(j).pitchPOST(targmed,POSTdays(i):POSTdays(i+1))-mean(PC24(j).pitchPRE(targmed,:)))/sigma,ws),'Color','k')
%         end
%     end
%     plot([-25 70],[0 0])
%     firstday=find(PC24(j).tvPOST<(PC24(j).tvPOST(1)+10));
%     avgadaptive(j,:)=coef*(mean(PC24(j).pitchPOST(targmed-32-200:targmed-32+200,firstday)')-mean(PC24(j).pitchPRE(targmed-32-200:targmed-32+200,:)'));
%     thr=0.5*coef;
%     plot([-25:5:70],[ones(1,20)*thr],'.')
%     xlim([-25 70])
%     ylim([-2 2])
% end





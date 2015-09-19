
%%%%%%
%%%%%%
% Example figure 
            load /bulbul1/LMAN_microstimulation/Data020912a.mat
        ExpIndex=4;
            figure;hold on;
            catchfrac=5;
            for i=36:45%length(experiment)
                window=[400:440];
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if endpoints(j)>startpoints(j)+catchfrac-1
                        inds=[startpoints(j):catchfrac:endpoints(j)];
                        plot(timedata(inds),FFdata(inds),'.','markersize',7,'Color','b')
                    end
                end
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if endpoints(j)>startpoints(j)+catchfrac-1
                        inds=[startpoints(j):catchfrac:endpoints(j)];           
                        plot(timedata(inds),FFdata(inds),'.','markersize',7,'Color','r')
                    end
                end
            end
            %%%%%%%%%%%%%%%%
            for i=1:length(IntExp(ExpIndex).FFmorningfitFBbaseline)
                t1=timing4(IntExp(ExpIndex).AllDaysFBbaseline(i).fv(1));
                t2=timing4(IntExp(ExpIndex).AllDaysFBbaseline(i).fv(end));
                plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBbaseline(i);IntExp(ExpIndex).FFeveningfitFBbaseline(i)],'r-','Linewidth',4)
            end

            for i=1:length(IntExp(ExpIndex).FFmorningfitFBlearning)
                t1=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(1));
                t2=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(end));
                plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBlearning(i);IntExp(ExpIndex).FFeveningfitFBlearning(i)],'r-','Linewidth',4)
            end
            %
            for i=1:length(IntExp(ExpIndex).FFmorningfitCTbaseline)
                t1=timing4(IntExp(ExpIndex).AllDaysCTbaseline(i).fv(1));
                t2=timing4(IntExp(ExpIndex).AllDaysCTbaseline(i).fv(end));
                plot([t1;t2],[IntExp(ExpIndex).FFmorningfitCTbaseline(i);IntExp(ExpIndex).FFeveningfitCTbaseline(i)],'g-','Linewidth',4)
            end

            for i=1:length(IntExp(ExpIndex).FFmorningfitCTlearning)
                t1=timing4(IntExp(ExpIndex).AllDaysCTlearning(i).fv(1));
                t2=timing4(IntExp(ExpIndex).AllDaysCTlearning(i).fv(end));
                plot([t1;t2],[IntExp(ExpIndex).FFmorningfitCTlearning(i);IntExp(ExpIndex).FFeveningfitCTlearning(i)],'g-','Linewidth',4)
            end
            plot([586 596],[mean(IntExp(ExpIndex).AllDaysCTbaseline(1).FFdata) mean(IntExp(ExpIndex).AllDaysCTbaseline(1).FFdata)],'g-','Linewidth',2)
            plot([586 596],[mean(IntExp(ExpIndex).AllDaysFBbaseline(1).FFdata) mean(IntExp(ExpIndex).AllDaysFBbaseline(1).FFdata)],'b-','Linewidth',2)

            plot([715 725],[mean(IntExp(ExpIndex).AllDaysCTlearning(5).FFdata) mean(IntExp(ExpIndex).AllDaysCTlearning(5).FFdata)],'g-','Linewidth',2)
            plot([715 725],[mean(IntExp(ExpIndex).AllDaysFBlearning(5).FFdata) mean(IntExp(ExpIndex).AllDaysFBlearning(5).FFdata)],'b-','Linewidth',2)
            learnday=0;
            learnnight=0;
            transferday=0;
            transfernight=0;
            for i=1:5
                learnday=learnday+(IntExp(ExpIndex).FFeveningfitCTlearning(i)-IntExp(ExpIndex).FFmorningfitCTlearning(i));
                 transferday=transferday+(IntExp(ExpIndex).FFeveningfitFBlearning(i)-IntExp(ExpIndex).FFmorningfitFBlearning(i));
            end
            for i=1:4
                learnnight=learnnight+(IntExp(ExpIndex).FFmorningfitCTlearning(i+1)-IntExp(ExpIndex).FFeveningfitCTlearning(i));
                transfernight=transfernight+(IntExp(ExpIndex).FFmorningfitFBlearning(i+1)-IntExp(ExpIndex).FFeveningfitFBlearning(i));  
            end

            plot([586 586],[mean(IntExp(ExpIndex).AllDaysCTbaseline(1).FFdata) mean(IntExp(ExpIndex).AllDaysCTbaseline(1).FFdata)+learnday],'r')
            plot([592 592],[mean(IntExp(ExpIndex).AllDaysCTbaseline(1).FFdata) mean(IntExp(ExpIndex).AllDaysCTbaseline(1).FFdata)+learnnight],'r')

            plot([715 715],[mean(IntExp(ExpIndex).AllDaysFBbaseline(1).FFdata) mean(IntExp(ExpIndex).AllDaysFBbaseline(1).FFdata)+transferday],'r')
            plot([725 725],[mean(IntExp(ExpIndex).AllDaysFBbaseline(1).FFdata) mean(IntExp(ExpIndex).AllDaysFBbaseline(1).FFdata)+transfernight],'r')


            xlim([581 727])                         
            ylim([2050 2450])
            
            
% Summary figures
            load /bulbul1/LMAN_microstimulation/Data020912a.mat
            
% if needed:            
edit stim020812.m


%%% All data - correlations

% In cases where there is much consolidation, was that due to big changes
% in the previous night or big changes over the course of the day?
                    % Stim trials
                            clear LearnST;
                            for i=1:6
                                % 1=all, 2=day, 3=night
                                LearnST{i}{1}=(IntExp(i).deltadayFBlearning(1:end-1)+IntExp(i).deltanightFBlearning(1:end));
                                LearnST{i}{2}=IntExp(i).deltadayFBlearning(1:end-1);
                                LearnST{i}{3}=IntExp(i).deltanightFBlearning(1:end);
                            end
                            clear allST dayST nightST
                            allST=([LearnST{2}{1} LearnST{3}{1} LearnST{4}{1}(1:4) LearnST{5}{1}(1:5)]);% LearnST{6}{1}]);
                            dayST=([LearnST{2}{2} LearnST{3}{2} LearnST{4}{2}(1:4) LearnST{5}{2}(1:5)]);% LearnST{6}{2}]);
                            nightST=([LearnST{2}{3} LearnST{3}{3} LearnST{4}{3}(1:4) LearnST{5}{3}(1:5)]);% LearnST{6}{3}]);
                    % Catch trials
                            clear LearnCT;
                            for i=1:6
                                % 1=all, 2=day, 3=night
                                LearnCT{i}{1}=(IntExp(i).deltadayCTlearning(1:end-1)+IntExp(i).deltanightCTlearning(1:end));
                                LearnCT{i}{2}=IntExp(i).deltadayCTlearning(1:end-1);
                                LearnCT{i}{3}=IntExp(i).deltanightCTlearning(1:end);
                            end
                            clear allCT dayCT nightCT
                            allCT=([LearnCT{2}{1} LearnCT{3}{1} LearnCT{4}{1}(1:4) LearnCT{5}{1}(1:5)]);% LearnCT{6}{1}]);
                            dayCT=([LearnCT{2}{2} LearnCT{3}{2} LearnCT{4}{2}(1:4) LearnCT{5}{2}(1:5)]);% LearnCT{6}{2}]);
                            nightCT=([LearnCT{2}{3} LearnCT{3}{3} LearnCT{4}{3}(1:4) LearnCT{5}{3}(1:5)]);% LearnCT{6}{3}]);
        Dirmatrix=[ones(1,9) -1*ones(1,4) ones(1,5)];

        
        
            figure;hold on;
            subplot(323);hold on;
            plot(dayST,allST,'ko')
%             p=polyfit(dayST,allST,1);
%             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
            xlim([-80 80]);ylim([-80 80])
            plot([-80 80],[0 0],'k');plot([0 0],[-80 80],'k')
            subplot(324);hold on;
            plot(nightST,allST,'ro')
%             p=polyfit(nightST,allST,1);
%             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
            xlim([-80 80]);ylim([-80 80])
             plot([-80 80],[0 0],'k');plot([0 0],[-80 80],'k')
           
            
%             figure;hold on;
%             subplot(121);hold on;
%             plot(dayST.*Dirmatrix,allST.*Dirmatrix,'ko')
% %             p=polyfit(dayST.*Dirmatrix,allST.*Dirmatrix,1);
% %             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
%             xlim([-80 80]);ylim([-80 80])
%             subplot(122);hold on;
%             plot(nightST.*Dirmatrix,allST.*Dirmatrix,'ro')
% %             p=polyfit(nightST.*Dirmatrix,allST.*Dirmatrix,1);
% %             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
%             xlim([-80 80]);ylim([-80 80])
            subplot(321);hold on;
            plot(dayCT,allCT,'ko')
%             p=polyfit(dayCT,allCT,1);
%             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
            xlim([-80 80]);ylim([-80 80])
            plot([-80 80],[0 0],'k');plot([0 0],[-80 80],'k')
            subplot(322);hold on;
            plot(nightCT,allCT,'ro')
%             p=polyfit(nightCT,allCT,1);
%             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
            xlim([-80 80]);ylim([-80 80])
             plot([-80 80],[0 0],'k');plot([0 0],[-80 80],'k')

            subplot(325);hold on;
            plot(dayCT-dayST,allCT-allST,'ko')
%             p=polyfit(dayCT,allCT,1);
%             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
            xlim([-80 80]);ylim([-80 80])
            plot([-80 80],[0 0],'k');plot([0 0],[-80 80],'k')
            subplot(326);hold on;
            plot(nightCT-nightST,allCT-allST,'ro')
%             p=polyfit(nightCT,allCT,1);
%             plot([-80:1:80],[-80:1:80]*p(1)+p(2))
            xlim([-80 80]);ylim([-80 80])
             plot([-80 80],[0 0],'k');plot([0 0],[-80 80],'k')
            
            
%%% During learning/transfer - bar plots

                            % Stim trials
                            lastday=[ 5 4 3 4 5 3];
                            Dirmatrix=[ones(1,7) -1*ones(1,4) ones(1,5)];
                            IntExp(1).dir=1;IntExp(2).dir=1;IntExp(3).dir=1;IntExp(4).dir=-1;IntExp(5).dir=1;IntExp(6).dir=-1;
                            clear LearnST;
                            for i=1:6
                                LearnST{i}{1}=IntExp(i).dir*(IntExp(i).deltadayFBlearning(1:lastday(i))+IntExp(i).deltanightFBlearning(1:lastday(i)));
                                LearnST{i}{2}=IntExp(i).dir*IntExp(i).deltadayFBlearning(1:lastday(i));
                                LearnST{i}{3}=IntExp(i).dir*IntExp(i).deltanightFBlearning(1:lastday(i));
                            end

                            allST=([LearnST{2}{1} LearnST{3}{1} LearnST{4}{1} LearnST{5}{1}]);% LearnST{6}{1}]);
                            dayST=([LearnST{2}{2} LearnST{3}{2} LearnST{4}{2} LearnST{5}{2}]);% LearnST{6}{2}]);
                            nightST=([LearnST{2}{3} LearnST{3}{3} LearnST{4}{3} LearnST{5}{3}]);% LearnST{6}{3}]);
                            % Catch trials
                            IntExp(1).dir=1;IntExp(2).dir=1;IntExp(3).dir=1;IntExp(4).dir=-1;
                            clear LearnCT;
                            for i=1:6
                                % 1=all, 2=day, 3=night
                                LearnCT{i}{1}=IntExp(i).dir*(IntExp(i).deltadayCTlearning(1:lastday(i))+IntExp(i).deltanightCTlearning(1:lastday(i)));
                                LearnCT{i}{2}=IntExp(i).dir*IntExp(i).deltadayCTlearning(1:lastday(i));
                                LearnCT{i}{3}=IntExp(i).dir*IntExp(i).deltanightCTlearning(1:lastday(i));
                            end

                            allCT=([LearnCT{2}{1} LearnCT{3}{1} LearnCT{4}{1} LearnCT{5}{1}]);% LearnCT{6}{1}]);
                            dayCT=([LearnCT{2}{2} LearnCT{3}{2} LearnCT{4}{2} LearnCT{5}{2}]);% LearnCT{6}{2}]);
                            nightCT=([LearnCT{2}{3} LearnCT{3}{3} LearnCT{4}{3} LearnCT{5}{3}]);% LearnCT{6}{3}]);
                            %
                figure;hold on;
                % total
                bar(1,mean(dayCT));bar(2,mean(nightCT))
                plot([1 1],[mean(dayCT)-std(dayCT)/sqrt(length(dayCT)) mean(dayCT)+std(dayCT)/sqrt(length(dayCT))],'-');
                plot([2 2],[mean(nightCT)-std(nightCT)/sqrt(length(nightCT)) mean(nightCT)+std(nightCT)/sqrt(length(nightCT))],'-');
                
                % motor pathway
                bar(4,mean(dayST));bar(5,mean(nightST))
                plot([4 4],[mean(dayST)-std(dayST)/sqrt(length(dayST)) mean(dayST)+std(dayST)/sqrt(length(dayST))],'-');
                plot([5 5],[mean(nightST)-std(nightST)/sqrt(length(nightST)) mean(nightST)+std(nightST)/sqrt(length(nightST))],'-');
                
                % AFP
                bar(7,mean(dayCT)-mean(dayST));bar(8,mean(nightCT)-mean(nightST))
                plot([7 7],[mean(dayCT-dayST)-std(dayCT-dayST)/sqrt(length(dayCT-dayST)) mean(dayCT-dayST)+std(dayCT-dayST)/sqrt(length(dayCT-dayST))],'-');
                plot([8 8],[mean(nightCT-nightST)-std(nightCT-nightST)/sqrt(length(nightCT-nightST)) mean(nightCT-nightST)+std(nightCT-nightST)/sqrt(length(nightCT-nightST))],'-');
                
%%% Just the first day


lastday=[ones(1,6)];
                            Dirmatrix=[ones(1,7) -1*ones(1,4) ones(1,5)];
                            IntExp(1).dir=1;IntExp(2).dir=1;IntExp(3).dir=1;IntExp(4).dir=-1;IntExp(5).dir=1;IntExp(6).dir=-1;
                            clear LearnST;
                            for i=1:6
                                LearnST{i}{1}=IntExp(i).dir*(IntExp(i).deltadayFBlearning(1:lastday(i))+IntExp(i).deltanightFBlearning(1:lastday(i)));
                                LearnST{i}{2}=IntExp(i).dir*IntExp(i).deltadayFBlearning(1:lastday(i));
                                LearnST{i}{3}=IntExp(i).dir*IntExp(i).deltanightFBlearning(1:lastday(i));
                            end

                            allST=([LearnST{2}{1} LearnST{3}{1} LearnST{4}{1} LearnST{5}{1}]);% LearnST{6}{1}]);
                            dayST=([LearnST{2}{2} LearnST{3}{2} LearnST{4}{2} LearnST{5}{2}]);% LearnST{6}{2}]);
                            nightST=([LearnST{2}{3} LearnST{3}{3} LearnST{4}{3} LearnST{5}{3}]);% LearnST{6}{3}]);
                            % Catch trials
                            IntExp(1).dir=1;IntExp(2).dir=1;IntExp(3).dir=1;IntExp(4).dir=-1;
                            clear LearnCT;
                            for i=1:6
                                % 1=all, 2=day, 3=night
                                LearnCT{i}{1}=IntExp(i).dir*(IntExp(i).deltadayCTlearning(1:lastday(i))+IntExp(i).deltanightCTlearning(1:lastday(i)));
                                LearnCT{i}{2}=IntExp(i).dir*IntExp(i).deltadayCTlearning(1:lastday(i));
                                LearnCT{i}{3}=IntExp(i).dir*IntExp(i).deltanightCTlearning(1:lastday(i));
                            end

                            allCT=([LearnCT{2}{1} LearnCT{3}{1} LearnCT{4}{1} LearnCT{5}{1}]);% LearnCT{6}{1}]);
                            dayCT=([LearnCT{2}{2} LearnCT{3}{2} LearnCT{4}{2} LearnCT{5}{2}]);% LearnCT{6}{2}]);
                            nightCT=([LearnCT{2}{3} LearnCT{3}{3} LearnCT{4}{3} LearnCT{5}{3}]);% LearnCT{6}{3}]);
                            %
                figure;hold on;
                % total
                bar(1,mean(dayCT));bar(2,mean(nightCT))
                plot([1 1],[mean(dayCT)-std(dayCT)/sqrt(length(dayCT)) mean(dayCT)+std(dayCT)/sqrt(length(dayCT))],'-');
                plot([2 2],[mean(nightCT)-std(nightCT)/sqrt(length(nightCT)) mean(nightCT)+std(nightCT)/sqrt(length(nightCT))],'-');
                
                % motor pathway
                bar(4,mean(dayST));bar(5,mean(nightST))
                plot([4 4],[mean(dayST)-std(dayST)/sqrt(length(dayST)) mean(dayST)+std(dayST)/sqrt(length(dayST))],'-');
                plot([5 5],[mean(nightST)-std(nightST)/sqrt(length(nightST)) mean(nightST)+std(nightST)/sqrt(length(nightST))],'-');
                
                % AFP
                bar(7,mean(dayCT)-mean(dayST));bar(8,mean(nightCT)-mean(nightST))
                plot([7 7],[mean(dayCT-dayST)-std(dayCT-dayST)/sqrt(length(dayCT-dayST)) mean(dayCT-dayST)+std(dayCT-dayST)/sqrt(length(dayCT-dayST))],'-');
                plot([8 8],[mean(nightCT-nightST)-std(nightCT-nightST)/sqrt(length(nightCT-nightST)) mean(nightCT-nightST)+std(nightCT-nightST)/sqrt(length(nightCT-nightST))],'-');

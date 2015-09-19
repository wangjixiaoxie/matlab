%%%% 
load /cardinal6/covertAnaly/Covert204.mat



%% ANALYSIS 1A - Is there learning?
                % 16 experiments --> 9 up, 8 down --> learning is much bigger for the up
                % experiments -- mean of 10Hz vs 50Hz

                % FFafter uses variable postday --- first and last song in day after exp
                        % postday --- rows are exps, column 1 is first point and column 2 is last point
                % For experiments 13 and 18:24, the postday is also the first time the bird
                % sings.  Mean learning in these experiments is 20Hz but n.s. by themselves

                % For all experiments that AP5 left and there was an opportunity for learning.
                % Exclude 3 because this was the pu56 one where too much AP5 caused song
                % degradation.  Exclude 8:9 because learning appeared to happen but
                % variability never recovered.  Exclude 11 because it's just a repeat of
                % 10 and thus is not an independent experiment (despite appearance of
                % covert learning).  

                Experiment=[Exps1 Exps2 Exps3];
                            %%%%%%%
                            %%% ONE EXTREME --- Baseline as all notes in pre file
                            %%%%%%%
                            clear Learning
                            ind=[1:2 4:7 10 12:24];
                            for j=1:length(ind)
                                i=ind(j);
                                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                                window=[Experiment(i).on:Experiment(i).off];
                                FFbefore=mean(mean(Experiment(i).pitchACpre(window,:)));
                                FFafter=mean(mean(Experiment(i).pitchACpost(window,postday(i,1):postday(i,2))));
                                Learning(i)=coef*(FFafter-FFbefore);
                            end
                            % p=<0.001 that Learning(ind) is different from zero - t-test
                            % mean and median are 29.3 and ..., respectively
                            % p=0.0127 for a sign test
                            
                            
                            %%%%%%%
                       % WHAT IF WE DON'T WAIT UNTIL THE NEXT DAY????
                            %%% ONE EXTREME --- Baseline as all notes in pre file
                            %%%%%%%
                            clear Learning
                            ind=[1:2 4:7 10 12:24];
                            for j=1:length(ind)
                                i=ind(j);
                                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                                window=[Experiment(i).on:Experiment(i).off];
                                FFbefore=mean(mean(Experiment(i).pitchACpre(window,:)));
                                Etime=Experiment(i).timeACpost(:,1:end);
                                [Esorted,sortedPost]=sort(Etime);
                                FFafter=mean(mean(Experiment(i).pitchACpost(window,sortedPost(1):postday(i,2))));
                                Learning(i)=coef*(FFafter-FFbefore);
                            end
                            % p=0.0014 that Learning(ind) is different from zero - t-test
                            % mean is 30.3, median is 24.2
                            % p=0.0127 for a sign test
                            % p=0.0197 if you just use the first 50 post
                            % p=0.0127 for a sign test
                            
                            
                            
                            %%%%%%%
                            %%% OTHER EXTREME ---- Baseline as last 20 notes in pre file
                            %%%%%%%
                            clear Learning
                            ind=[1:2 4:7 10 12:24];
                            for j=1:length(ind)
                                i=ind(j);
                                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                                window=[Experiment(i).on:Experiment(i).off];
                                FFbefore=mean(mean(Experiment(i).pitchACpre(window,end-20:end)));
                                FFafter=mean(mean(Experiment(i).pitchACpost(window,postday(i,1):postday(i,2))));
                                Learning(i)=coef*(FFafter-FFbefore);
                            end
                            % p=0.0064
                            % mean and median are 34.2 and 27.6, respectively
                            % p=0.0127 for a sign test

                            %%%%%%
                            %%% EXTREME OF BOTH --- Baseline as last 20 notes in pre file, Post as
                            %%% first 20 notes in post file --- 20 notes gives SEM of about 10Hz
                            %%%%%%
                            clear Learning
                            ind=[1:2 4:7 10:24];
                            for j=1:length(ind)
                                i=ind(j);
                                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                                window=[Experiment(i).on:Experiment(i).off];
                                FFbefore=mean(mean(Experiment(i).pitchACpre(window,end-20:end)));
                                FFafter=mean(mean(Experiment(i).pitchACpost(window,postday(i,1):postday(i,1)+20)));
                                Learning(i)=coef*(FFafter-FFbefore);
                            end
                            [h,p]=ttest(Learning(ind))
                            % p=0.0212
                            % Still significant, but barely. If last 50 and first 50, then it's 0.05 
                            % mean and median are 29.5 and 28.9, respectively
                            % p=0.21 for a sign test --- n.s.
                                % p=0.02 for a sign test with all pre data and first 100
                                % post points, p=0.07 for a sign test with last 20 pre and
                                % first 100 post
%%% ANALYSIS 1B - Is learning syllable-specific? 
            % This controls for effect of AP5, effect of non-specific WN anxiety, etc.
                        clear ControlNotes
                        postday2=postday;
                        postday2(15,:)=[3 24];
                        postday2(19,:)=[1 93]
                        postday2(21,:)=[1 20];
                        postday2(22,:)=[1 13];  
                        postday2(23,:)=[1 103];                        
                        ind2=[1:2 7 10 12:24];
                        for j=1:length(ind2)
                            i=ind2(j);
                            coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                            window=[Experiment(i).onCTL:Experiment(i).offCTL];
                            FFbefore(i)=mean(mean(Experiment(i).pitchACpreCTL(window,:)));
                            FFafter(i)=mean(mean(Experiment(i).pitchACpostCTL(window,postday2(i,1):postday2(i,2))));
                           %% for ones where postday is not the same as for
                           %% the Targeted note
                           FFapv(i)=coef*(mean(mean(Experiment(i).pitchAPVCTL(window,round(end/2):end)))-FFbefore(i));
                            LearningCTL(i)=coef*(FFafter(i)-FFbefore(i));
                        end
                        % p=0.7007 -- from t-test -- definitely syllable specific
                
%%% ANALYSIS 1C - Does AP5 really block learning? 
            %%% CHECK FOR FFchange with AP5 on ---
            for i=1:length(ind)
                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                window=[Experiment(i).on:Experiment(i).off];
                % All data from the beginning of the next day onwards
                Etime=Experiment(i).timeAPVwn;
                [x,sortedAPVwn]=sort(Etime);
                Etime2=Experiment(i).timeAPV;
                [x,sortedAPV]=sort(Etime2);

                Epitch=Experiment(i).pitchAPVwn;
                Epitch2=Experiment(i).pitchAPV;
                endAPV(i)=coef*(mean(mean(Epitch(window,sortedAPVwn(end/5:end))))-mean(mean(Experiment(i).pitchACpre(window,:))));
                beginAPV(i)=coef*(mean(mean(Epitch2(window,sortedAPV(end-20:end))))-mean(mean(Experiment(i).pitchACpre(window,:))));

            end
            AP5n=endAPV-beginAPV;
            % mean is 11.76, p=0.26 with test 
            
       %%% The significant difference between these values and the learning in a paired t-test
          % is somewhat sensitive to how "endAPV" is defined --- i.e. last
          % twenty percent of syllables
                figure;hold on
                plot([1;2],[beginAPV(ind);endAPV(ind)])
                plot(1,beginAPV(ind),'V','Markersize',10,'Color','k')
                plot(2,endAPV(ind),'V','Markersize',10,'Color','k')    
                plot([0 3],[0 0],'k')
                xlim([0.7 2.3])   
                plot([0.8 1.2],[mean(beginAPV(ind)) mean(beginAPV(ind))],'-','Color','k','Linewidth',3)
                plot([1.8 2.2],[mean(endAPV(ind)) mean(endAPV(ind))],'-','Color','k','Linewidth',3)  
                ylim([-100 100])
        %%% DEFINE APV end as the last hour with significant amount of song
                    for i=ind
                            coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                            window=[Experiment(i).on:Experiment(i).off];
                            % All data from the beginning of the next day onwards
                            Etime=Experiment(i).timeAPVwn;
                            [xWN,sortedAPVwn]=sort(Etime);
                            Etime2=Experiment(i).timeAPV;
                            [x,sortedAPV]=sort(Etime2);
                            APVwnvals=find(xWN>max(xWN(end-10))-1);
                            Epitch=Experiment(i).pitchAPVwn;
                            Epitch2=Experiment(i).pitchAPV;           
                            endAPV(i)=coef*(mean(mean(Epitch(window,APVwnvals)))-mean(mean(Experiment(i).pitchACpre(window,:))));
                            beginAPV(i)=coef*(mean(mean(Epitch2(window,end-20:end)))-mean(mean(Experiment(i).pitchACpre(window,:))));
                    end
                figure;hold on
               % plot([1;2],[beginAPV(ind);endAPV(ind)])
                plot(1,beginAPV(ind),'V','Markersize',10,'Color','k')
                %plot(2,endAPV(ind),'V','Markersize',10,'Color','k')    
                plot([0 3],[0 0],'k')
                xlim([0.7 1.3])   
                plot([0.8 1.2],[mean(beginAPV(ind)) mean(beginAPV(ind))],'-','Color','k','Linewidth',3)
                %plot([1.8 2.2],[mean(endAPV(ind)) mean(endAPV(ind))],'-','Color','k','Linewidth',3)  
                ylim([-100 100])
        
        
%%% ANALYSIS 1D - Positive controls                    
                % See below - analysis 5
% FIGURE 1 - Controls vs. Experimental
            figure;plot(1,Learning(find(Learning~=0)),'V','Color','k','Markersize',10);
            hold on;plot(3,LearningCTL(find(abs(LearningCTL)>0)),'V','Color','b','Markersize',10)     
            hold on;plot(0,AP5n(find(abs(AP5n)>0)),'V','Color','r','Markersize',10)     
            plot(2,finalpointPC,'V','Color','c','Markersize',10)
            plot([-0.5 3.5],[0 0],'k')
            xlim([-0.5 3.5])
            plot([0.7 1.3],[mean(Learning(find(abs(Learning)>0))) mean(Learning(find(abs(Learning)>0)))],'Linewidth',3)
            plot([2.7 3.3],[mean(LearningCTL(find(abs(LearningCTL)>0))) mean(LearningCTL(find(abs(LearningCTL)>0)))],'Linewidth',3)
            plot([-0.3 0.3],[mean(AP5n(find(abs(AP5n)>0))) mean(AP5n(find(abs(AP5n)>0)))],'Linewidth',3)
            plot([1.7 2.3],[mean(finalpointPC) mean(finalpointPC)],'Linewidth',3)
            ylim([-100 100])                
%%% ANALYSIS 2 - aesthetically pleasing examples
    % FIGURE 2
    % bk13bk12
            window=[300:360];
            figure;hold on;
            plot(Experiment(14).timeACpre,(mean(Experiment(14).pitchACpre(window,:))-mean(mean(Experiment(14).pitchACpre(window,:)))),'.','Markersize',10)
            plot(Experiment(14).timeAPV(11:end),(mean(Experiment(14).pitchAPV(window,11:end))-mean(mean(Experiment(14).pitchACpre(window,:)))),'.','Color','g','Markersize',10)
            plot(Experiment(14).timeAPVwn,(mean(Experiment(14).pitchAPVwn(window,:))-mean(mean(Experiment(14).pitchACpre(window,:)))),'.','Color','r','Markersize',10)
            plot(Experiment(14).timeACpost(250:382),(mean(Experiment(14).pitchACpost(window,250:382))-mean(mean(Experiment(14).pitchACpre(window,:)))),'.','Color','k','Markersize',10)
            avgescapes=mean(mean(pitchE(window,:)))-mean(mean(Experiment(14).pitchACpre(window,:)));
            plot([8165 8180],[mean(mean(Experiment(14).pitchACpost(window,250:382)))-mean(mean(Experiment(14).pitchACpre(window,:))) mean(mean(Experiment(14).pitchACpost(window,250:382)))-mean(mean(Experiment(14).pitchACpre(window,:)))] ,'Color','k')
            plot([8148 8155],[avgescapes avgescapes])
            xlim([8120 8180])
            ylim([-120 120])
            plot([8120 8180],[0 0],'k')
    % bk75bk62 --- Experiment #12-13
             figure;hold on;
                avg927escapes=mean(mean(pitchY(250:400,find(std(pitchY(400:500,:))<50))'))-mean(mean(Experiment(12).pitchACpre(250:400,36:end)));
                avg1002escapes=mean(mean(pitch1002escape(250:400,:)))-mean(mean(Experiment(13).pitchACpre(250:400,:)));  
                % 36:end because PRE is 24 hrs - control for circadian!
             subplot(221);hold on;plot(Experiment(12).timeACpre(36:111),mean(Experiment(12).pitchACpre(250:400,36:111))-mean(mean(Experiment(12).pitchACpre(250:400,36:end))),'.','Markersize',10)
                plot(Experiment(12).timeACpre(112:end),mean(Experiment(12).pitchACpre(250:400,112:end))-mean(mean(Experiment(12).pitchACpre(250:400,36:end))),'.','Markersize',10)
                % 33:end because APV off before dark
                hold on;plot(Experiment(12).timeACpost(1:269),mean(Experiment(12).pitchACpost(250:400,1:269))-mean(mean(Experiment(12).pitchACpre(250:400,36:end))),'.','Markersize',10,'Color','k')
                hold on;plot(-217.5+Experiment(12).timeAPV,mean(Experiment(12).pitchAPV(250:400,:))-mean(mean(Experiment(12).pitchACpre(250:400,36:end))),'.','Markersize',10,'Color','g')
                hold on;plot(-217.5+Experiment(12).timeAPVwn,mean(Experiment(12).pitchAPVwn(250:400,:))-mean(mean(Experiment(12).pitchACpre(250:400,36:end))),'.','Markersize',10,'Color','r') 
                hold on;plot(runningaverage(Experiment(12).timeACpost(1:269),50),runningaverage(mean(Experiment(12).pitchACpost(250:400,1:269))-mean(mean(Experiment(12).pitchACpre(250:400,36:end))),50),'Color','r','Linewidth',3)
                plot([7117 7163],[0 0],'k')
                xlim([7117 7163]);ylim([-200 200])
                plot([7149 7157],[avg927escapes avg927escapes])
                % 114:end
             subplot(222);hold on;plot(Experiment(13).timeACpre(1:54),mean(Experiment(13).pitchACpre(250:400,1:54))-mean(mean(Experiment(13).pitchACpre(250:400,:))),'.','Markersize',10)
                hold on;plot(Experiment(13).timeACpre(55:end),mean(Experiment(13).pitchACpre(250:400,55:end))-mean(mean(Experiment(13).pitchACpre(250:400,:))),'.','Markersize',10)
                hold on;plot(Experiment(13).timeACpost(1:end),mean(Experiment(13).pitchACpost(250:400,1:end))-mean(mean(Experiment(13).pitchACpre(250:400,:))),'.','Markersize',10,'Color','k')
                hold on;plot(Experiment(13).timeAPV,mean(Experiment(13).pitchAPV(250:400,:))-mean(mean(Experiment(13).pitchACpre(250:400,:))),'.','Markersize',10,'Color','g')
                hold on;plot(Experiment(13).timeAPVwn,mean(Experiment(13).pitchAPVwn(250:400,:))-mean(mean(Experiment(13).pitchACpre(250:400,:))),'.','Markersize',10,'Color','r')
                hold on;plot(runningaverage(Experiment(13).timeACpost(1:end),50),runningaverage(mean(Experiment(13).pitchACpost(250:400,1:end))-mean(mean(Experiment(13).pitchACpre(250:400,:))),50),'Color','r','Linewidth',3)
                plot([7242 7283],[0 0],'k')
                xlim([7242 7283]);ylim([-200 200])
                plot([7268 7276],[avg1002escapes avg1002escapes])
             subplot(223);hold on;plot(-217+Experiment(12).timeACpreCTL(1:76),mean(Experiment(12).pitchACpreCTL(220:340,1:76))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),'.','Markersize',10)
                plot(-217+Experiment(12).timeACpreCTL(77:end),mean(Experiment(12).pitchACpreCTL(220:340,77:end))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),'.','Markersize',10)
                % 33:end because APVCTL off before dark
                hold on;plot(-217+Experiment(12).timeACpostCTL(1:end),mean(Experiment(12).pitchACpostCTL(220:340,1:end))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),'.','Markersize',10,'Color','k')
                hold on;plot(-217+Experiment(12).timeAPVCTL,mean(Experiment(12).pitchAPVCTL(220:340,:))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),'.','Markersize',10,'Color','g')
                hold on;plot(-217+Experiment(12).timeAPVwnCTL,mean(Experiment(12).pitchAPVwnCTL(220:340,:))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),'.','Markersize',10,'Color','r') 
                hold on;plot(-217+runningaverage(Experiment(12).timeACpostCTL(1:end),50),runningaverage(mean(Experiment(12).pitchACpostCTL(230:340,1:end))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),50),'Color','r','Linewidth',3)
                plot([7117 7163],[0 0],'k')
                xlim([7117 7163]);ylim([-200 200])
                % 114:end
            subplot(224);plot(672+Experiment(13).timeACpreCTL,mean(Experiment(13).pitchACpreCTL(220:340,:))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),'.','Markersize',10)
                hold on;plot(672+Experiment(13).timeACpostCTL(1:end),mean(Experiment(13).pitchACpostCTL(220:340,1:end))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),'.','Markersize',10,'Color','k')
                plot(672+Experiment(13).timeAPVCTL,mean(Experiment(13).pitchAPVCTL(220:340,:))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),'.','Markersize',10,'Color','g')
                plot(672+Experiment(13).timeAPVwnCTL,mean(Experiment(13).pitchAPVwnCTL(220:340,:))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),'.','Markersize',10,'Color','r')
                hold on;plot(672+runningaverage(Experiment(13).timeACpostCTL(1:end),50),runningaverage(mean(Experiment(13).pitchACpostCTL(230:340,1:end))-mean(mean(Experiment(13).pitchACpreCTL(220:340,1:end))),50),'Color','r','Linewidth',3)
                plot([7242 7283],[0 0],'k')
                xlim([7242 7283]);ylim([-200 200])
                        %     % Experiments 19 and 20 display bidirectionality in a nice fashion despite
                        %     % 20 not actually being significant learning based on the criteria
                        %     % above.
                        %    figure;hold on;
                        %             % 36:end because PRE is 24 hrs - control for circadian!
                        %  subplot(211);hold on;plot(Experiment(12).timeACpre(36:end),mean(Experiment(12).pitchACpre(320:380,36:end))-mean(mean(Experiment(12).pitchACpre(320:380,36:end))),'*')
                        %             % 33:end because APV off before dark
                        %     hold on;plot(Experiment(12).timeACpost(33:end),mean(Experiment(12).pitchACpost(320:380,33:end))-mean(mean(Experiment(12).pitchACpre(320:380,36:end))),'*','Color','k')
                        %     hold on;plot(Experiment(12).timeAPV,mean(Experiment(12).pitchAPV(320:380,:))-mean(mean(Experiment(12).pitchACpre(320:380,36:end))),'*','Color','g')
                        %     hold on;plot(Experiment(12).timeAPVwn,mean(Experiment(12).pitchAPVwn(320:380,:))-mean(mean(Experiment(12).pitchACpre(320:380,36:end))),'*','Color','r') 
                        %     hold on;plot(Experiment(13).timeACpre,mean(Experiment(13).pitchACpre(320:380,:))-mean(mean(Experiment(13).pitchACpre(320:380,:))),'*')
                        %     hold on;plot(Experiment(13).timeACpost(114:end),mean(Experiment(13).pitchACpost(320:380,114:end))-mean(mean(Experiment(13).pitchACpre(320:380,:))),'*','Color','k')
                        %     hold on;plot(Experiment(13).timeAPV,mean(Experiment(13).pitchAPV(320:380,:))-mean(mean(Experiment(13).pitchACpre(320:380,:))),'*','Color','g')
                        %     hold on;plot(Experiment(13).timeAPVwn,mean(Experiment(13).pitchAPVwn(320:380,:))-mean(mean(Experiment(13).pitchACpre(320:380,:))),'*','Color','r')
                        %     plot([7100 7300],[0 0],'k')
                        %     ylim([-220 220])
                        %     xlim([7117 7288])
                        %  subplot(212);hold on;plot(Experiment(12).timeACpre(36:end),mean(Experiment(12).pitchACpreCTL(1250:1350,36:end))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,36:end))),'*')
                        %             % 33:end because APVCTL off before dark
                        %     hold on;plot(Experiment(12).timeACpostCTL(33:end),mean(Experiment(12).pitchACpostCTL(1250:1350,33:end))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,36:end))),'*','Color','k')
                        %     hold on;plot(Experiment(12).timeAPV,mean(Experiment(12).pitchAPVCTL(1250:1350,:))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,36:end))),'*','Color','g')
                        %     hold on;plot(Experiment(12).timeAPVwn,mean(Experiment(12).pitchAPVwnCTL(1250:1350,:))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,36:end))),'*','Color','r') 
                        %     hold on;plot(Experiment(13).timeACpreCTL,mean(Experiment(13).pitchACpreCTL(1250:1350,:))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*')
                        %     hold on;plot(Experiment(13).timeACpostCTL(114:end),mean(Experiment(13).pitchACpostCTL(1250:1350,114:end))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*','Color','k')
                        %     hold on;plot(Experiment(13).timeAPV,mean(Experiment(13).pitchAPVCTL(1250:1350,:))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*','Color','g')
                        %     hold on;plot(Experiment(13).timeAPVwn,mean(Experiment(13).pitchAPVwnCTL(1250:1350,:))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*','Color','r')
                        %     plot([7100 7300],[0 0],'k')
                        %     ylim([-220 220])
                        %     xlim([7117 7288])    
                        %     
                        %        figure;hold on;
                        %             % 36:end because PRE is 24 hrs - control for circadian!
                        %     subplot(211);hold on;plot(runningaverage(Experiment(12).timeACpre(36:111),20),runningaverage(mean(Experiment(12).pitchACpre(250:350,36:111))-mean(mean(Experiment(12).pitchACpre(250:350,36:end))),20),'*')
                        %     plot(runningaverage(Experiment(12).timeACpre(112:end),20),runningaverage(mean(Experiment(12).pitchACpre(250:350,112:end))-mean(mean(Experiment(12).pitchACpre(250:350,36:end))),20),'*')
                        %     % 33:end because APV off before dark
                        %     hold on;plot(runningaverage(Experiment(12).timeACpost(33:269),20),runningaverage(mean(Experiment(12).pitchACpost(250:350,33:269))-mean(mean(Experiment(12).pitchACpre(250:350,36:end))),20),'*','Color','k')
                        %     hold on;plot(-217.5+runningaverage(Experiment(12).timeAPV,20),runningaverage(mean(Experiment(12).pitchAPV(250:350,:))-mean(mean(Experiment(12).pitchACpre(250:350,36:end))),20),'*','Color','g')
                        %     hold on;plot(-217.5+runningaverage(Experiment(12).timeAPVwn,20),runningaverage(mean(Experiment(12).pitchAPVwn(250:350,:))-mean(mean(Experiment(12).pitchACpre(250:350,36:end))),20),'*','Color','r') 
                        %     hold on;plot(runningaverage(Experiment(13).timeACpre(1:54),20),runningaverage(mean(Experiment(13).pitchACpre(250:350,1:54))-mean(mean(Experiment(13).pitchACpre(250:350,:))),20),'*')
                        %     hold on;plot(runningaverage(Experiment(13).timeACpre(55:end),20),runningaverage(mean(Experiment(13).pitchACpre(250:350,55:end))-mean(mean(Experiment(13).pitchACpre(250:350,:))),20),'*')
                        %     hold on;plot(runningaverage(Experiment(13).timeACpost(164:end),20),runningaverage(mean(Experiment(13).pitchACpost(250:350,164:end))-mean(mean(Experiment(13).pitchACpre(250:350,:))),20),'*','Color','k')
                        %     hold on;plot(runningaverage(Experiment(13).timeAPV,20),runningaverage(mean(Experiment(13).pitchAPV(250:350,:))-mean(mean(Experiment(13).pitchACpre(250:350,:))),20),'*','Color','g')
                        %     hold on;plot(runningaverage(Experiment(13).timeAPVwn,20),runningaverage(mean(Experiment(13).pitchAPVwn(250:350,:))-mean(mean(Experiment(13).pitchACpre(250:350,:))),20),'*','Color','r')
                        %     plot([7100 7300],[0 0],'k')
                        %     ylim([-100 100])
                        % 
                        %  subplot(212);hold on;plot(-217+runningaverage(Experiment(12).timeACpreCTL(1:76),20),runningaverage(mean(Experiment(12).pitchACpreCTL(220:340,1:76))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),20),'*')
                        %  plot(-217+runningaverage(Experiment(12).timeACpreCTL(77:end),20),runningaverage(mean(Experiment(12).pitchACpreCTL(220:340,77:end))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),20),'*')
                        %  % 33:end because APVCTL off before dark
                        %     hold on;plot(-217+runningaverage(Experiment(12).timeACpostCTL(33:end),20),runningaverage(mean(Experiment(12).pitchACpostCTL(220:340,33:end))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),20),'*','Color','k')
                        %     hold on;plot(-217+runningaverage(Experiment(12).timeAPVCTL,20),runningaverage(mean(Experiment(12).pitchAPVCTL(220:340,:))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),20),'*','Color','g')
                        %     hold on;plot(-217+runningaverage(Experiment(12).timeAPVwnCTL,20),runningaverage(mean(Experiment(12).pitchAPVwnCTL(220:340,:))-mean(mean(Experiment(12).pitchACpreCTL(220:340,1:end))),20),'*','Color','r') 
                        %     hold on;plot(672+runningaverage(Experiment(13).timeACpreCTL,20),runningaverage(mean(Experiment(13).pitchACpreCTL(220:340,:))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),20),'*')
                        %     hold on;plot(672+runningaverage(Experiment(13).timeACpostCTL(61:end),20),runningaverage(mean(Experiment(13).pitchACpostCTL(220:340,61:end))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),20),'*','Color','k')
                        %     plot(672+runningaverage(Experiment(13).timeAPVCTL,20),runningaverage(mean(Experiment(13).pitchAPVCTL(220:340,:))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),20),'*','Color','g')
                        %     plot(672+runningaverage(Experiment(13).timeAPVwnCTL,20),runningaverage(mean(Experiment(13).pitchAPVwnCTL(220:340,:))-mean(mean(Experiment(13).pitchACpreCTL(220:340,:))),20),'*','Color','r')
                        %     plot([7100 7300],[0 0],'k')
                        %     ylim([-100 100])
                        % 
%%% ANALYSIS 3 - Does AP5 block learning by disrupting LMAN-->RA?
%%% CHECKS FOR VARIABILITY REDUCTION and RECOVERY - JUSTIFIES "next day"

               % APV --- ROBUST to defining APV as last 20 or 30 or last 25% in APV batch
                 % Regardless which ones you exclude, results remain around p=0.003
                 % So...don't exclude any!
                 % On average, var reduction is around 25-30%

               % POST --- is defined as first 73 to allow #22, which doesn't have much data
                 % On average, post is around 97% of pre

              % More notes about robust-ness:
              % (VarAPV-VarPre)./VarPre -- All are > 17% except #16,19
              %                         -- Conclusions remain similarly significant when you
              %                             exclude 4,16,19 - all with low inas -- b/c
              %                             although 4 was +60, 16 was -30Hz(maladaptive)
              % (VarPost-VarPre)./VarPre
        figure;hold on;
        % ANALYSIS 3A - Targeted notes, targeted window, CV     
             for j=1:length(ind)
                 i=ind(j);
                 window3=[Experiment(i).on:Experiment(i).off];
                % This is just because some parts are noisy - doesn't
                % affect outcome.  Var is reduced everywhere
                 if i==16;window3=[250:350];end
                 if i==4;window3=[180:230];end
                 if i==7;window3=[250:350];end
                 if i==10;window3=[300:400];end
                 a=mean(std(Experiment(i).pitchACpre(window3,:)'));
                 b=mean(mean(Experiment(i).pitchACpre(window3,:)));
                 VarPre(i)=a/b;
                 a=mean(std(Experiment(i).pitchAPV(window3,end-20:end)'));
                 b=mean(mean(Experiment(i).pitchAPV(window3,end-20:end)));                 
                 VarAPV(i)=a/b;
                 a=mean(std(Experiment(i).pitchACpost(window3,postday(i,1):postday(i,1)+100)'));
                 b=mean(mean(Experiment(i).pitchACpost(window3,postday(i,1):postday(i,1)+100)));
                 VarPost(i)=a/b;
                 if i==19
                    a=mean(std(Experiment(i).pitchACpreCTL(300:400,:)'));
                    b=mean(mean(Experiment(i).pitchACpreCTL(300:400,:)));
                    VarPre(i)=a/b;
                    a=mean(std(pitchAPV(300:400,round(end-20):end)'))
                    b=mean(mean(pitchAPV(300:400,round(end-20):end)));
                    VarAPV(i)=a/b;
                    a=mean(std(Experiment(i).pitchACpostCTL(300:400,postday(i,1):postday(i,1)+100)'));
                    b=mean(mean(Experiment(i).pitchACpostCTL(300:400,postday(i,1):postday(i,1)+100)));
                    VarPost(i)=a/b;
                 end
             end
           subplot(131);hold on;
           plot([1;2;3],[VarPre(ind);VarAPV(ind);VarPost(ind)])  
           plot([0.8 1.2],[mean(VarPre(ind)) mean(VarPre(ind))],'Color','k','Linewidth',3)
           plot([1.8 2.2],[mean(VarAPV(ind)) mean(VarAPV(ind))],'Color','k','Linewidth',3)
           plot([2.8 3.2],[mean(VarPost(ind)) mean(VarPost(ind))],'Color','k','Linewidth',3)
           xlim([0.6 3.4])
           ylim([0 0.035])
        % ANALYSIS 3B - Targeted notes, targeted window, Absolute values     
             for j=1:length(ind)
                 i=ind(j);
                 window=[Experiment(i).on:Experiment(i).off];
                 if i==16;window=[250:350];end
                 if i==4;window=[180:230];end
                 if i==7;window=[250:350];end
                 if i==10;window=[300:400];end                 
                 VarPre(i)=mean(std(Experiment(i).pitchACpre(window,:)'));
                 VarAPV(i)=mean(std(Experiment(i).pitchAPV(window,end-20:end)'));
                 VarPost(i)=mean(std(Experiment(i).pitchACpost(window,postday(i,1):postday(i,1)+50)'));
                 if i==19 %%% All other notes show variability constancy. Choose CTL notes
                    VarPre(i)=mean(std(Experiment(i).pitchACpreCTL(300:400,:)'));
                    VarAPV(i)=mean(std(pitchAPV(300:400,round(end-20):end)'));
                    VarPost(i)=mean(std(Experiment(i).pitchACpostCTL(300:400,postday(i,1):postday(i,1)+50)'));
                 end

             end
           subplot(132);hold on;
           plot([1;2;3],[VarPre(ind);VarAPV(ind);VarPost(ind)])  
           plot([0.8 1.2],[mean(VarPre(ind)) mean(VarPre(ind))],'Color','k','Linewidth',3)
           plot([1.8 2.2],[mean(VarAPV(ind)) mean(VarAPV(ind))],'Color','k','Linewidth',3)
           plot([2.8 3.2],[mean(VarPost(ind)) mean(VarPost(ind))],'Color','k','Linewidth',3)
           xlim([0.6 3.4])
           ylim([0 130])

        % ANALYSIS 3C - Targeted notes, targeted window, Normalized
           subplot(133);hold on;
           plot([1;2;3],[VarPre(ind)./VarPre(ind);VarAPV(ind)./VarPre(ind);VarPost(ind)./VarPre(ind)])  
           plot([0.8 1.2],[1 1],'Color','k','Linewidth',3)
           plot([1.8 2.2],[mean(VarAPV(ind))/mean(VarPre(ind)) mean(VarAPV(ind))/mean(VarPre(ind))],'Color','k','Linewidth',3)
           plot([2.8 3.2],[mean(VarPost(ind))/mean(VarPre(ind)) mean(VarPost(ind))/mean(VarPre(ind))],'Color','k','Linewidth',3)
           xlim([0.6 3.4])
           ylim([0 2])

        % I still need to quantify APV data for most of the control notes
        % ANALYSIS 3D - Control notes, control window, CV
        
        
   %%% What if we don't wait until overnight     
        figure;hold on;
        % ANALYSIS 3E - Targeted notes, targeted window, CV     
             for j=1:length(ind)
                 i=ind(j);
                 window=[Experiment(i).on:Experiment(i).off];
                % There are some noisy parts of these notes 
                    % This is NOT statistical chicanery.
                 if i==16;window=[250:350];end
                 if i==4;window=[180:230];end
                 if i==7;window=[250:350];end
                 if i==10;window=[300:400];end
                 Etime=Experiment(i).timeACpost(:,1:end);
                [Esorted,sortedPost]=sort(Etime);
                 VarPre(i)=mean(std(Experiment(i).pitchACpre(window,:)'))/mean(mean(Experiment(i).pitchACpre(window,:)));
                 VarAPV(i)=mean(std(Experiment(i).pitchAPV(window,end-20:end)'))/mean(mean(Experiment(i).pitchAPV(window,end-20:end)));
                 VarPost(i)=mean(std(Experiment(i).pitchACpost(window,sortedPost(1:50))'))/mean(mean(Experiment(i).pitchACpost(window,sortedPost(1:50))));
                 if i==19
                    VarPre(i)=mean(std(Experiment(i).pitchACpreCTL(300:400,:)'))/mean(mean(Experiment(i).pitchACpreCTL(300:400,:)));
                    VarAPV(i)=mean(std(pitchAPV(300:400,round(end-20):end)'))/mean(mean(pitchAPV(300:400,round(end-20):end)));
                    VarPost(i)=mean(std(Experiment(i).pitchACpostCTL(300:400,1:50)'))/mean(mean(Experiment(i).pitchACpostCTL(300:400,1:50)));
                 end
             
             end
             
           subplot(131);hold on;
           plot([1;2;3],[VarPre(ind);VarAPV(ind);VarPost(ind)])  
           plot([0.8 1.2],[mean(VarPre(ind)) mean(VarPre(ind))],'Color','k','Linewidth',3)
           plot([1.8 2.2],[mean(VarAPV(ind)) mean(VarAPV(ind))],'Color','k','Linewidth',3)
           plot([2.8 3.2],[mean(VarPost(ind)) mean(VarPost(ind))],'Color','k','Linewidth',3)
           xlim([0.6 3.4])
           ylim([0 0.035])

%%% ANALYSIS 4 - Decay/Recovery/Extinction/Forgetting analysis       
% Figure 4A - Total hours (includes sleep)

        figure;hold on
        subplot(221);hold on
        for i=ind
%             colorchoice(i,:)=[rand rand rand];
            coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
            window=[Experiment(i).on:Experiment(i).off];
            % All data from the beginning of the next day onwards
            Etime=Experiment(i).timeACpost(:,postday(i,1):end);
            [x,sortedPost]=sort(Etime);
            Epitch=Experiment(i).pitchACpost(:,postday(i,1):end);
            plot(runningaverage(Etime(sortedPost)-Etime(sortedPost(1)),20),...
                coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(Experiment(i).pitchACpre(window,:))),20)),'-','Color',colorchoice(i,:))
            timevals(i).data=Etime(sortedPost)-Etime(sortedPost(1));
            pitchvals(i).data=coef*(mean(Epitch(window,sortedPost))-mean(mean(Experiment(i).pitchACpre(window,:))));
        end
        % CALCULATE MEAN DURING EACH HOUR WHERE n>2
        for i=1:79
            allvals=[];
            for j=ind
                hourvals=find(timevals(j).data>i-1 & timevals(j).data<i);
                if ~isempty(pitchvals(j).data(hourvals))
                    allvals=[allvals mean(pitchvals(j).data(hourvals))];
                end
            end
            meanFF(i)=mean(allvals);
        end
        plot(meanFF,'Color','b','Linewidth',4)
        plot([0 100],[0 0],'Color','k')
        xlim([0 100])
        ylim([-150 150])
% Figure 4B - Hours awake
clear diff;
    subplot(222);hold on
    for i=ind
            coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
            window=[Experiment(i).on:Experiment(i).off];
            Etime=Experiment(i).timeACpost(:,postday(i,1):end);
            [Esorted,sortedPost]=sort(Etime);
          % Finds indices of gaps longer than 8hrs
            nights=find(diff(Esorted)>8);
          % Subtracts gaps
            for j=1:length(nights)
                gaplength=Esorted(nights(j)+1)-Esorted(nights(j));
                Esorted(nights(j)+1:end)=Esorted(nights(j)+1:end)-gaplength;
            end
          % plots points with appropriate x-axis
            Epitch=Experiment(i).pitchACpost(:,postday(i,1):end);
            plot(runningaverage(Esorted-Esorted(1),20),...
                coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(Experiment(i).pitchACpre(window,:))),20)),'Color',colorchoice(i,:))
            timevalsNOGAP(i).data=Esorted-Esorted(1);
    end
            % CALCULATE MEAN DURING EACH HOUR WHERE n>2
        for i=1:40
            allvals=[];
            for j=ind
                hourvals=find(timevalsNOGAP(j).data>i-1 & timevalsNOGAP(j).data<i);
                if ~isempty(pitchvals(j).data(hourvals))
                    allvals=[allvals mean(pitchvals(j).data(hourvals))];
                end
            end
            meanFFnogap(i)=mean(allvals);
            stdFFnogap(i)=std(allvals);
            nums(i)=length(allvals)
        end
        plot(meanFFnogap,'Color','b','Linewidth',4)
        plot([1:40;1:40],[meanFFnogap-stdFFnogap./sqrt(nums-1);meanFFnogap+stdFFnogap./sqrt(nums-1)],'Color','b','Linewidth',4)
        plot([0 100],[0 0],'Color','k')
        xlim([0 45])
        ylim([-150 150])
%%% FOR C and D, we also look at time before the first night
        % Figure 4C - Total hours (includes sleep)
                subplot(223);hold on
                for i=ind
                    coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                    window=[Experiment(i).on:Experiment(i).off];
                    % All data from the beginning of the next day onwards
                    Etime=Experiment(i).timeACpost(:,1:end);
                    [x,sortedPost]=sort(Etime);
                    Epitch=Experiment(i).pitchACpost(:,1:end);
                    plot(runningaverage(Etime(sortedPost)-Etime(sortedPost(1)),20),...
                        coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(Experiment(i).pitchACpre(window,:))),20)),'-','Color',colorchoice(i,:))
                    timevals(i).data=Etime(sortedPost)-Etime(sortedPost(1));
                    pitchvals(i).data=coef*(mean(Epitch(window,sortedPost))-mean(mean(Experiment(i).pitchACpre(window,:))));
                end
                % CALCULATE MEAN DURING EACH HOUR WHERE n>2
                for i=1:79
                    allvals=[];
                    for j=ind
                        hourvals=find(timevals(j).data>i-1 & timevals(j).data<i);
                        if ~isempty(pitchvals(j).data(hourvals))
                            allvals=[allvals mean(pitchvals(j).data(hourvals))];
                        end
                    end
                    meanFF(i)=mean(allvals);
                end
                plot(meanFF,'Color','b','Linewidth',4)
                plot([0 100],[0 0],'Color','k')
                xlim([0 100])
                ylim([-150 150])
        % Figure 4D - Hours awake
            subplot(224);hold on
            for i=ind
                    coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                    window=[Experiment(i).on:Experiment(i).off];
                    Etime=Experiment(i).timeACpost(:,1:end);
                    [Esorted,sortedPost]=sort(Etime);
                  % Finds indices of gaps longer than 8hrs
                    nights=find(diff(Esorted)>8);
                  % Subtracts gaps
                    for j=1:length(nights)
                        gaplength=Esorted(nights(j)+1)-Esorted(nights(j));
                        Esorted(nights(j)+1:end)=Esorted(nights(j)+1:end)-gaplength;
                    end
                  % plots points with appropriate x-axis
                    Epitch=Experiment(i).pitchACpost(:,1:end);
                    plot(runningaverage(Esorted-Esorted(1),20),...
                        coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(Experiment(i).pitchACpre(window,:))),20)),'Color',colorchoice(i,:))
                    timevalsNOGAP(i).data=Esorted-Esorted(1);
            end
                    % CALCULATE MEAN DURING EACH HOUR WHERE n>2
                for i=1:40
                    allvals=[];
                    for j=ind
                        hourvals=find(timevalsNOGAP(j).data>i-1 & timevalsNOGAP(j).data<i);
                        if ~isempty(pitchvals(j).data(hourvals))
                            allvals=[allvals mean(pitchvals(j).data(hourvals))];
                        end
                    end
                    meanFFnogap(i)=mean(allvals);
                    stdFFnogap(i)=std(allvals);
                    nums(i)=length(allvals)
                end
                plot(meanFFnogap,'Color','b','Linewidth',4)
                plot([1:40;1:40],[meanFFnogap-stdFFnogap./sqrt(nums-1);meanFFnogap+stdFFnogap./sqrt(nums-1)],'Color','b','Linewidth',4)
                plot([0 100],[0 0],'Color','k')
                xlim([0 45])
                ylim([-150 150])


%%% ANALYSIS 5 - Red nucleus krupa/thompson figure
            figure;hold on;
            subplot(312);hold on;
            for i=ind
                window=[Experiment(i).on:Experiment(i).off];
                EtimeAPV=Experiment(i).timeAPVwn;
                [Esorted,sortedAPV]=sort(EtimeAPV);
                mpitchpre=mean(mean(Experiment(i).pitchACpre(window,:)));
                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                onethird(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(end/4)-5:round(end/4)+5))))...
                    -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                twothirds(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(2*end/4)-5:round(2*end/4)+5))))...
                    -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                finalpoint(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(end-10:end))))...
                    -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                post1(i)=(coef*(mean(mean(Experiment(i).pitchACpost(window,1:20)))...
                    -mean(mean(Experiment(i).pitchACpre(window,:)))))/mpitchpre;
                post2(i)=(coef*(mean(mean(Experiment(i).pitchACpost(window,1:postday(i,2)-1)))...
                    -mean(mean(Experiment(i).pitchACpre(window,:)))))/mpitchpre;
                LearningNorm(ind)=Learning(ind)/mpitchpre;

            end
            denom1=sqrt(length(ind));

            errorbar([1,2,3,4,6,7],100*[0,mean(onethird(ind)),mean(twothirds(ind)),mean(finalpoint(ind)),mean(post1(ind)),mean(post2(ind))],...
                100*[0,std(onethird(ind))/denom1,std(twothirds(ind))/denom1,std(finalpoint(ind))/denom1,std(post1(ind))/denom1,std(post2(ind))/denom1],'Color','r')

            for i=1:length(ExperimentPC)
                coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
                mpitchpre=mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)));
                onethirdPC(i)=(coef*(mean(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,round(end/4)-5:round(end/4)+5)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                twothirdsPC(i)=(coef*(mean(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,round(2*end/4)-5:round(2*end/4)+5)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                finalpointPC(i)=(coef*(mean(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,end-10:end)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                postPC1(i)=(coef*(mean(mean(ExperimentPC(i).pitchPost(ExperimentPC(i).time,1:20)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                postPC2(i)=(coef*(mean(mean(ExperimentPC(i).pitchPost(ExperimentPC(i).time,:)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                postPCff(i).data=(coef*((mean(ExperimentPC(i).pitchPost(ExperimentPC(i).time-10:ExperimentPC(i).time+10,:)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                postPCtime(i).data=(ExperimentPC(i).timePost)-ExperimentPC(i).timePost(1);

            end

%             for i=1:14
%                 for j=2:5
%                     thishour(j,i)=mean(postPCff(i).data(find(floor(postPCtime(i).data)==j)));
%                 end
%             end
%             figure;hold on;
            
            
            denom=sqrt(length(ExperimentPC));
            errorbar([1,2,3,4,6,7],100*[0,mean(onethirdPC),mean(twothirdsPC),mean(finalpointPC),mean(postPC1),mean(postPC2)],...
                100*[0,std(onethirdPC)/denom,std(twothirdsPC)/denom,std(finalpointPC)/denom,std(postPC1)/denom,std(postPC2)/denom])
            plot([0.5 7.5],[0 0],'k')
            xlim([0.5 7.5])
            ylim([-1 2])


            %%%% 
            % WHAT IF WE JUST LOOKED AT THE ONES WITHOUT SIGNIFICANT AP5 LEARNING
            % (<15Hz)
            acceptable=[1 4 7 9 10 12 13 14 15 17 18 19];
            denom1=sqrt(length((ind(acceptable))));
            subplot(313);hold on;
            errorbar([1,2,3,4,6,7],100*[0,mean(onethird(ind(acceptable))),mean(twothirds(ind(acceptable))),mean(finalpoint(ind(acceptable))),mean(post1(ind(acceptable))),mean(post2(ind(acceptable)))],...
                100*[0,std(onethird(ind(acceptable)))/denom1,std(twothirds(ind(acceptable)))/denom1,std(finalpoint(ind(acceptable)))/denom1,std(post1(ind(acceptable)))/denom1,std(post2(ind(acceptable)))/denom1],'Color','r')
            errorbar([1,2,3,4,6,7],100*[0,mean(onethirdPC),mean(twothirdsPC),mean(finalpointPC),mean(postPC1),mean(postPC2)],...
                100*[0,std(onethirdPC)/denom,std(twothirdsPC)/denom,std(finalpointPC)/denom,std(postPC1)/denom,std(postPC2)/denom])
            plot([0.5 7.5],[0 0],'k')
            xlim([0.5 7.5])
            ylim([-1 2])
subplot(311);hold on;
plot(100*finalpoint(ind),100*LearningNorm(ind),'.','Markersize',15)
plot(100*finalpoint(ind(acceptable)),100*LearningNorm(ind(acceptable)),'.','Markersize',15,'Color','r')
plot(100*[-0.035 0.035],[0 0],'k')
plot([0 0],100*[-0.035 0.035],'k')
xlim(100*[-0.035 0.035])
ylim(100*[-0.035 0.035])

            %%%%%%%
            % WHAT IF WE JUST LOOKED AT THE ONES WITH A TON OF VARIABILITY
            % REDUCTION???
            figure;hold on;
            subplot(311);hold on;
            plot(1-(VarAPV(ind)./VarPre(ind)),finalpoint(ind)*100,'.','Markersize',15)
            %acceptable=[1 4 7 9 10 13 14 19];
            %plot(1-(VarAPV(ind(acceptable))./VarPre(ind(acceptable))),finalpoint(ind(acceptable))*100,'.','Markersize',15,'Color','r')
            plot([-0.5 1],[0 0],'k')
            plot([0 0],[-2.2 2.2],'k')
            %plot([0.3 0.3],[-2.2 2.2],'k')
            ylim([-2.2 2.2])
            subplot(312);hold on;
            acceptable=[1 2 4 5 7 8 9 10 13 14 19 20];
            denom1=sqrt(length((ind(acceptable))));
            errorbar([1,2,3,4,6,7],100*[0,mean(onethird(ind(acceptable))),mean(twothirds(ind(acceptable))),mean(finalpoint(ind(acceptable))),mean(post1(ind(acceptable))),mean(post2(ind(acceptable)))],...
                100*[0,std(onethird(ind(acceptable)))/denom1,std(twothirds(ind(acceptable)))/denom1,std(finalpoint(ind(acceptable)))/denom1,std(post1(ind(acceptable)))/denom1,std(post2(ind(acceptable)))/denom1],'Color','r')
            errorbar([1,2,3,4,6,7],100*[0,mean(onethirdPC),mean(twothirdsPC),mean(finalpointPC),mean(postPC1),mean(postPC2)],...
                100*[0,std(onethirdPC)/denom,std(twothirdsPC)/denom,std(finalpointPC)/denom,std(postPC1)/denom,std(postPC2)/denom])
            plot([0.5 7.5],[0 0],'k')
            xlim([0.5 7.5])
            ylim([-1 2])

            %%%% 
            % WHAT IF WE JUST LOOKED AT THE ONES WITH BOTH VAR REDUCTION
            % AND BLOCK OF EXPRESSION???
            subplot(313);hold on
            acceptable=[1 4 7 9 10 13 14 19];
            denom1=sqrt(length((ind(acceptable))));
            errorbar([1,2,3,4,6,7],100*[0,mean(onethird(ind(acceptable))),mean(twothirds(ind(acceptable))),mean(finalpoint(ind(acceptable))),mean(post1(ind(acceptable))),mean(post2(ind(acceptable)))],...
                100*[0,std(onethird(ind(acceptable)))/denom1,std(twothirds(ind(acceptable)))/denom1,std(finalpoint(ind(acceptable)))/denom1,std(post1(ind(acceptable)))/denom1,std(post2(ind(acceptable)))/denom1],'Color','r')
            errorbar([1,2,3,4,6,7],100*[0,mean(onethirdPC),mean(twothirdsPC),mean(finalpointPC),mean(postPC1),mean(postPC2)],...
                100*[0,std(onethirdPC)/denom,std(twothirdsPC)/denom,std(finalpointPC)/denom,std(postPC1)/denom,std(postPC2)/denom])
            plot([0.5 7.5],[0 0],'k')
            xlim([0.5 7.5])
            ylim([-1 2]) 
            


            



%% ANALYSIS 6: Example #2
    figure;hold on;
    window=[Experiment(i).on:Experiment(i).off];
    coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
    plot(runningaverage(Experiment(2).timeACpre(211:463),20),runningaverage(mean(Experiment(2).pitchACpre(window,211:463))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10)
    plot(runningaverage(Experiment(2).timeACpre(464:end),20),runningaverage(mean(Experiment(2).pitchACpre(window,464:end))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10)
    plot(runningaverage(Experiment(2).timeAPV,20),runningaverage(mean(Experiment(2).pitchAPV(window,:))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10,'Color','g')
    plot(runningaverage(Experiment(2).timeAPVwn(1:132),20),runningaverage(mean(Experiment(2).pitchAPVwn(window,1:132))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10,'Color','r')
    plot(runningaverage(Experiment(2).timeAPVwn(140:end),20),runningaverage(mean(Experiment(2).pitchAPVwn(window,140:end))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10,'Color','r')
    plot(runningaverage(Experiment(2).timeACpost(1:postday(2,1)),20),runningaverage(mean(Experiment(2).pitchACpost(window,1:postday(2,1)))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10,'Color','k')
    plot(runningaverage(Experiment(2).timeACpost(postday(2,1)+1:postday(2,2)),20),runningaverage(mean(Experiment(2).pitchACpost(window,postday(2,1)+1:postday(2,2)))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10,'Color','k')
    plot(runningaverage(Experiment(2).timeACpost(postday(2,2)+1:1054),20),runningaverage(mean(Experiment(2).pitchACpost(window,postday(2,2)+1:1054))-mean(mean(Experiment(2).pitchACpre(window,211:end))),20),'.','Markersize',10,'Color','k')

    
%%%% ANALYSIS 7 - Localization within targeted note?

            change=[];
            figure;hold on
            for i=ind
                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
               changed=(coef*(median((Experiment(i).pitchACpost(nbounds(i,1):nbounds(i,2),1:postday(i,2))'))-mean(Experiment(i).pitchACpre(nbounds(i,1):nbounds(i,2),:)')))
               plot(changed)
               normchanged(i).data=changed;
            end
            centered=zeros(20,1000)
            for i=1:20
                j=ind(i);
                mnwin=median(Experiment(j).TargetingWN-32);
                diff=round(mnwin-nbounds(j,1));
                if diff>length(normchanged(j).data)
                    diff=length(normchanged(j).data);
                end
                firsthalf=normchanged(j).data(1:diff);
                secondhalf=normchanged(j).data(diff:end);
                centered(i,500-diff:499)=firsthalf;
                centered(i,500:499+length(secondhalf))=secondhalf;
            end
            for k=400:600
                mntrajectory(k-399)=mean(centered(find(centered(:,k)~=0),k));
            end
            figure;hold on;
            plot(mntrajectory)
            xlim([0 200])
            ylim([0 30])

   %%% What about Positive Controls? Is lack of localization due to AP5 block or short-term learning? 

            figure;hold on
            for i=[1:14]
               mtargPC=median(ExperimentPC(i).time)-32;
                coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
               changed=(coef*(median((ExperimentPC(i).pitchPost(nboundsPC(i,1):nboundsPC(i,2),1:20)'))-mean(ExperimentPC(i).pitchPre(nboundsPC(i,1):nboundsPC(i,2),:)')))
               if i==2 | i==7 | i==8
                    changed=(coef*(median((ExperimentPC(i).pitchWN(nboundsPC(i,1):nboundsPC(i,2),1:20)'))-mean(ExperimentPC(i).pitchPre(nboundsPC(i,1):nboundsPC(i,2),:)')))
               end
               plot(changed)
               normchangedPC(i).data=changed;
            end
            centeredPC=zeros(14,1000)
            for j=[1:14]
                mnwin=round(median(ExperimentPC(j).time-32));
                if j==9
                    mnwin=200;
                end
                diff=round(mnwin-nboundsPC(j,1));
                if diff>length(normchangedPC(j).data)
                    diff=length(normchangedPC(j).data);
                end
                firsthalf=normchangedPC(j).data(1:diff);
                secondhalf=normchangedPC(j).data(diff:end);
                centeredPC(j,500-diff:499)=firsthalf;
                centeredPC(j,500:499+length(secondhalf))=secondhalf;
            end
            clear mntrajectoryPC
            for k=400:600
                mntrajectoryPC(k-399)=mean(centeredPC(find(centeredPC(:,k)~=0),k));
            end
            figure;hold on;
            plot(mntrajectoryPC)
            xlim([20 180])
            ylim([0 35])
%%%%%%%%%%
mean(mean(Experiment(10).pitchACpre(Experiment(10).on:Experiment(10).off,end-50:end)))-mean(mean(Experiment(10).pitchAPV(Experiment(10).on:Experiment(10).off,round(end/2):end)))
mean(mean(Experiment(11).pitchACpre(Experiment(10).on:Experiment(10).off,end-50:end)))-mean(mean(Experiment(11).pitchAPV(Experiment(10).on:Experiment(10).off,round(end/2):end)))
    
    %%% OLD ANALYSIS
ind([1 4 6 8 10 11 12 14 15 16])



figure
FFrecov=[260 450 1 50 50 50 0 0 0 50 0 50 140 1 1 40 1 1 1 1 1 1 1 1];
for i=[1:2 4:6 10 12 13 15:16 18:24]
    target=round(median(Experiment(i).on:Experiment(i).off));
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    [x,sorted]=sort(Experiment(i).timeACpre);
    t1=[1/length(sorted):1/length(sorted):1];
    baseline=median(Experiment(i).pitchACpre(target,:));
    if i==16
        baseline=median(Experiment(i).pitchACpre(target,:));
    end
    hold on;plot(runningaverage(t1,round(length(t1)/20)),coef*(runningaverage(Experiment(i).pitchACpre(target,sorted),round(length(t1)/20))-baseline),'-')
    ratime(i,1).data=runningaverage(t1,round(length(t1)/20));
    raFF(i,1).data=coef*(runningaverage(Experiment(i).pitchACpre(target,sorted),round(length(t1)/20))-baseline);
    [x,sorted2]=sort(Experiment(i).timeAPV);
    t2=[1+1/length(sorted2):1/length(sorted2):2];
    hold on;plot(runningaverage(t2,round(length(t2)/20)),coef*(runningaverage(Experiment(i).pitchAPV(target,sorted2),round(length(t2)/20))-baseline),'-','Color','g')
    ratime(i,2).data=runningaverage(t2,round(length(t2)/20));
    raFF(i,2).data=coef*(runningaverage(Experiment(i).pitchAPV(target,sorted2),round(length(t2)/20))-baseline);
    [x,sorted3]=sort(Experiment(i).timeAPVwn);
    t3=[2+1/length(sorted3):1/length(sorted3):3];
    hold on;plot(runningaverage(t3,round(length(t3)/20)),coef*(runningaverage(Experiment(i).pitchAPVwn(target,sorted3),round(length(t3)/20))-baseline),'-','Color','r')
    ratime(i,3).data=runningaverage(t3,round(length(t3)/20));
    raFF(i,3).data=coef*(runningaverage(Experiment(i).pitchAPVwn(target,sorted3),round(length(t3)/20))-baseline);
    [x,sorted4]=sort(Experiment(i).timeACpost(FFrecov(i):end));
    t4=[3+1/length(sorted4):1/length(sorted4):4];
    time4(i).data=x;
    hold on;plot(runningaverage(t4,round(length(t4)/20)),coef*(runningaverage(Experiment(i).pitchACpost(target,sorted4),round(length(t4)/20))-baseline),'-','Color','k')
    ratime(i,4).data=runningaverage(t4,round(length(t4)/20));
    raFF(i,4).data=coef*(runningaverage(Experiment(i).pitchACpost(target,sorted4),round(length(t4)/20))-baseline);
end


ind=[1:2 4:6 10 12 13:16 18:24];
for i=[1:2 4:6 10 12 13:16 18:24]
    for j=1:20
        width=length(raFF(i,1).data)/20;
        mnFF1(i,j)=mean(raFF(i,1).data(1+(j-1)*width:j*width));
        width=length(raFF(i,2).data)/20;
        mnFF2(i,j)=mean(raFF(i,2).data(1+(j-1)*width:j*width));
        width=length(raFF(i,4).data)/20;
        mnFF4(i,j)=mean(raFF(i,4).data(1+(j-1)*width:j*width));
        width=length(raFF(i,3).data)/20;
        mnFF3(i,j)=mean(raFF(i,3).data(1+(j-1)*width:j*width));
    end
end
figure;hold on;
plot([1:20],mean(mnFF1(ind,:)))
plot([21:40],mean(mnFF2(ind,:)),'g')
plot([41:60],mean(mnFF3(ind,:)),'r')
plot([61:80],mean(mnFF4(ind,:)),'k')



ind=[1:2 4:6 10 12:17 18:24]
for i=ind
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    targ=Experiment(i).on;
    baseline=mean(Experiment(i).pitchACpre(targ,find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-5)));
    [x,sortedAPV]=sort(Experiment(i).timeAPV);
    APVbegin(i)=coef*(mean(Experiment(i).pitchAPV(targ,sortedAPV(end-15:end)))-baseline);
    APVend(i)=coef*(mean(Experiment(i).pitchAPVwn(targ,end-15:end))-baseline);
    for j=1:70
        times1=find(Experiment(i).timeACpost>min(Experiment(i).timeACpost)+j-1);
        times2=find(Experiment(i).timeACpost(times1)<min(Experiment(i).timeACpost)+j);
        place=[];
        place=coef*(Experiment(i).pitchACpost(targ,times1(times2))-baseline);
        if ~isempty(place)
            ACpost(i,j)=mean(place);
        else
            ACpost(i,j)=0;
        end
    end
end
for i=1:70
    mnACpost(i)=mean(ACpost(find(ACpost(:,i)~=0),i));
    sdACpost(i)=std(ACpost(find(ACpost(:,i)~=0),i));
end

figure;plot(mnACpost+sdACpost/sqrt(12))
hold on;plot(mnACpost-sdACpost/sqrt(12))
hold on;plot([-10 0],[mean(APVbegin(ind));mean(APVend(ind))])
hold on;plot([-10 0],[mean(APVbegin(ind))+std(APVbegin(ind))/sqrt(12);mean(APVend(ind))+std(APVend(ind))/sqrt(12)])
hold on;plot([-10 0],[mean(APVbegin(ind))-std(APVbegin(ind))/sqrt(12);mean(APVend(ind))-std(APVend(ind))/sqrt(12)])




figure;plot([APVbegin;APVend])

figure;hold on
for i=ind %[1:2 4:6 10 12:13 15:16 18:24]
if isequal(Experiment(i).DIR,'up');
coef=1;
else
coef=-1;
end
hold on;plot(runningaverage(Experiment(i).timeACpost,15)-min(Experiment(i).timeACpost),coef*(runningaverage(Experiment(i).pitchACpost(Experiment(i).on,:),15)-mean(Experiment(i).pitchACpre(Experiment(i).on,:))),'k')
end



%%%%
% with sleep
%%%%
ind=[1:2 4:6 10 12:13 15:16 18:24]
for i=ind
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    targ=Experiment(i).on;
    baseline=mean(Experiment(i).pitchACpre(targ,find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-5)));
    [x,sortedAPV]=sort(Experiment(i).timeAPV);
    APVbegin(i)=coef*(mean(Experiment(i).pitchAPV(targ,sortedAPV(end-15:end)))-baseline);
    APVend(i)=coef*(mean(Experiment(i).pitchAPVwn(targ,end-15:end))-baseline);
    % find first night of sleep
    [x,sortACpost]=sort(Experiment(i).timeACpost);
    clear gaps
    for k=2:length(Experiment(i).timeACpost)
        gaps(k-1)=Experiment(i).timeACpost(sortACpost(k))-Experiment(i).timeACpost(sortACpost(k-1));
    end
    nights=find(gaps>8);
    nighttime1=round((Experiment(i).timeACpost(sortACpost(nights(1)-1))-min(Experiment(i).timeACpost)));
    g(i)=nighttime1;
    for j=1:nighttime1
        times1=find(Experiment(i).timeACpost>min(Experiment(i).timeACpost)+j-1);
        times2=find(Experiment(i).timeACpost(times1)<min(Experiment(i).timeACpost)+j);
        place=[];
        place=coef*(Experiment(i).pitchACpost(targ,times1(times2))-baseline);
        if ~isempty(place)
            ACpost(i,j)=mean(place);
        else
            ACpost(i,j)=0;
        end
    end
    for j=nighttime1:70
        times1=find(Experiment(i).timeACpost>min(Experiment(i).timeACpost)+j-1);
        times2=find(Experiment(i).timeACpost(times1)<min(Experiment(i).timeACpost)+j);
        place=[];
        place=coef*(Experiment(i).pitchACpost(targ,times1(times2))-baseline);
        if ~isempty(place)
            ACpost2(i,j)=mean(place);
        else
            ACpost2(i,j)=0;
        end

    end
    
end
for i=1:70
    mnACpost2(i)=mean(ACpost2(find(ACpost2(:,i)~=0),i));
    sdACpost2(i)=std(ACpost2(find(ACpost2(:,i)~=0),i));
end

figure;plot(mnACpost+sdACpost/sqrt(12))
hold on;plot(mnACpost-sdACpost/sqrt(12))
hold on;plot([-10 0],[mean(APVbegin(ind));mean(APVend(ind))])
hold on;plot([-10 0],[mean(APVbegin(ind))+std(APVbegin(ind))/sqrt(12);mean(APVend(ind))+std(APVend(ind))/sqrt(12)])
hold on;plot([-10 0],[mean(APVbegin(ind))-std(APVbegin(ind))/sqrt(12);mean(APVend(ind))-std(APVend(ind))/sqrt(12)])

%%%%%%%%%%
%%%%%%%%%%
%% 12.9.09
%%%%%%%%%%
% Figure 1 - Adaptive shifts (entire next day post - pvs 50 songs pre)
figure;plot(1,shiftvals,'v','Color','b','MarkerSize',10)
hold on;plot([0 2],[0 0],'k')
hold on;plot([0.9 1.1],[mean(shiftvals) mean(shiftvals)],'k')
xlim([0.8 1.2])
hold on;plot(1,shiftvals([2 4 7 10 13]),'v','Color','r','MarkerSize',10)
% Figure 1prime - Adaptive shifts 2 
%Same as above but with pre being all data in the pre file - better and just as significant
shiftvals2=shiftvals;
shiftvals2(1:10)=shiftvals2(1:10)-differ([1:2 4:6 10 12:13 15:16]);
figure;plot(1,shiftvals2,'v','Color','b','MarkerSize',10)
hold on;plot([0 2],[0 0],'k')
hold on;plot([0.9 1.1],[mean(shiftvals2) mean(shiftvals2)],'k')
xlim([0.8 1.2])
hold on;plot(1,shiftvals2([2 4 7 10 13]),'v','Color','r','MarkerSize',10)
%%% What is the change with AP5 on? - subtract end of AP5 post (last 20%) - from end of AP5 w/o WN (last 20%)
    for i=[1:2 4:6 10 12:13 15:16]
    AP5n(i)=mean(raFF(i,3).data(end-length(raFF(i,3).data)/5:end))-mean(raFF(i,2).data(end-length(raFF(i,2).data)/5:end));
    end
    AP5n(17)=mean(pitch1201apvwn(240,end-size(pitch1201apvwn,2)/5:end))-mean(pitch1201apv(240,end-size(pitch1201apv,2)/5:end));
    AP5n(18)=mean(pitch1203apvwnC(260,end-size(pitch1203apvwnC,2)/5:end))-mean(pitch1203apvC(260,end-size(pitch1203apvC,2)/5:end));
    AP5n(19)=-1*(mean(pitch1206apvwnC(260,end-size(pitch1206apvwnC,2)/5:end))-mean(pitch1206apvC(260,end-size(pitch1206apvC,2)/5:end)));
    AP5n=AP5n([1:2 4:6 10 12:13 15:19]);
hold on;plot(0.5,AP5n,'v','Color','b','MarkerSize',10)
xlim([0.3 1.2])
plot([0.4 0.6],[mean(AP5n) mean(AP5n)],'k')
plot(0.5,AP5n([2 4 7 10 13]),'v','Color','r','MarkerSize',10)
plot([0.5;1],[AP5n;shiftvals2],'k')


%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%% 12.15.09 %%%

%%% Control notes - AP5 and WN but no learning signal %%%
% Experiments 1:2 10 12:13 16 use the next note in the repeat as the ctl.
% Experiment 15 uses the long low stack (whereas the short high stack is
    % targeted).
% Experiments 18:20 use the low stack or high stack (whichever one is not
    % targeted).
ind=[1:2 10 12:13 15 16 18:20]
for i=ind
    if isequal(Experiment(i).DIR,'up');
        coef=1;
    else
        coef=-1;
    end
    targ=Experiment(i).onCTL;
    baseline=mean(Experiment(i).pitchACpreCTL(targ,find(Experiment(i).timeACpreCTL>Experiment(i).timeACpreCTL(end)-5)));
    % find first night of sleep
        [x,sortACpost]=sort(Experiment(i).timeACpostCTL);
        clear gaps
        gaps(1)=Experiment(i).timeACpostCTL(sortACpost(1))-Experiment(i).timeAPVwn(end);
        for k=2:length(Experiment(i).timeACpostCTL)
            gaps(k)=Experiment(i).timeACpostCTL(sortACpost(k))-Experiment(i).timeACpostCTL(sortACpost(k-1));
        end
        nights=find(gaps>8);
        if length(nights)==1
            indexed=sortACpost(nights(1):end);
        else
            indexed=sortACpost(nights(1):nights(2));
        end
        diffCTL(i)=coef*(mean(Experiment(i).pitchACpostCTL(targ,indexed))-baseline);
end
for i=1:70
    mnACpost2(i)=mean(ACpost2(find(ACpost2(:,i)~=0),i));
    sdACpost2(i)=std(ACpost2(find(ACpost2(:,i)~=0),i));
end

%%%%%%%%%%
%%%%%%%%%%
for i=[1:2 4:6 10 12:13 15:16 18:24]
    

   figure;plot(Experiment(12).timeACpre,mean(Experiment(12).pitchACpreCTL(320:380,:))-mean(mean(Experiment(12).pitchACpreCTL(320:380,:))),'*')
            % 33:end because APVCTL off before dark
    hold on;plot(Experiment(12).timeACpostCTL(33:end),mean(Experiment(12).pitchACpostCTL(320:380,33:end))-mean(mean(Experiment(12).pitchACpreCTL(320:380,:))),'*','Color','k')
    hold on;plot(Experiment(12).timeAPVCTL,mean(Experiment(12).pitchAPVCTL(320:380,:))-mean(mean(Experiment(12).pitchACpreCTL(320:380,:))),'*','Color','g')
    hold on;plot(Experiment(12).timeAPVwnCTL,mean(Experiment(12).pitchAPVwnCTL(320:380,:))-mean(mean(Experiment(12).pitchACpreCTL(320:380,:))),'*','Color','r') 
    hold on;plot(Experiment(13).timeACpreCTL,mean(Experiment(13).pitchACpreCTL(320:380,:))-mean(mean(Experiment(13).pitchACpreCTL(320:380,:))),'*')
    hold on;plot(Experiment(13).timeACpostCTL(114:end),mean(Experiment(13).pitchACpostCTL(320:380,114:end))-mean(mean(Experiment(13).pitchACpreCTL(320:380,:))),'*','Color','k')
    hold on;plot(Experiment(13).timeAPVCTL,mean(Experiment(13).pitchAPVCTL(320:380,:))-mean(mean(Experiment(13).pitchACpreCTL(320:380,:))),'*','Color','g')
    hold on;plot(Experiment(13).timeAPVwnCTL,mean(Experiment(13).pitchAPVwnCTL(320:380,:))-mean(mean(Experiment(13).pitchACpreCTL(320:380,:))),'*','Color','r')
   figure;plot(Experiment(12).timeACpre,mean(Experiment(12).pitchACpreCTL(1250:1350,:))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,:))),'*')
            % 33:end because APVCTL off before dark
    hold on;plot(Experiment(12).timeACpostCTL(33:end),mean(Experiment(12).pitchACpostCTL(1250:1350,33:end))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,:))),'*','Color','k')
    hold on;plot(Experiment(12).timeAPV,mean(Experiment(12).pitchAPVCTL(1250:1350,:))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,:))),'*','Color','g')
    hold on;plot(Experiment(12).timeAPVwn,mean(Experiment(12).pitchAPVwnCTL(1250:1350,:))-mean(mean(Experiment(12).pitchACpreCTL(1250:1350,:))),'*','Color','r') 
    hold on;plot(Experiment(13).timeACpreCTL,mean(Experiment(13).pitchACpreCTL(1250:1350,:))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*')
    hold on;plot(Experiment(13).timeACpostCTL(114:end),mean(Experiment(13).pitchACpostCTL(1250:1350,114:end))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*','Color','k')
    hold on;plot(Experiment(13).timeAPV,mean(Experiment(13).pitchAPVCTL(1250:1350,:))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*','Color','g')
    hold on;plot(Experiment(13).timeAPVwn,mean(Experiment(13).pitchAPVwnCTL(1250:1350,:))-mean(mean(Experiment(13).pitchACpreCTL(1250:1350,:))),'*','Color','r')

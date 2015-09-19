load /cardinal5/Covert040810.mat



% 811 - figs
figure;hold on;
            plot([0.1 0.1],[pcCIbase],'b')
            plot([3.0 3.0],[pcCIfinal],'b')
            plot([3.4 3.4],[pcCIpost1],'b')
%             plot([4.9 4.9],[pcCIpost2],'b')
            plot([0.1 3.0 3.4],[0 mean(finalpointPC) mean(postPC1)],'b')
            plot([0 0],[CIbase],'r')
%             plot([1 1],[CIfirst],'r')
%             plot([2 2],[CIsecond],'r')
            %plot([3 3],[CIthird],'r')
            plot([3 3],[CIfinal],'r')
            plot([3.3 3.3],[CIpost1],'r')
%             plot([5 5],[CIpost2],'r')
            plot([0 3 3.3],[0 mean(finalpoint(ind)) mean(post1(ind)) ],'r')
            plot([-1 6],[0 0],'k')
figure;hold on;
bar([mean(VarPre(ind)) 0 0  mean(VarAPV(ind)) 0 0  mean(VarPost(ind))])
plot([1 1],[mean(VarPre(ind))-std(VarPre(ind))/length(ind) mean(VarPre(ind))+std(VarPre(ind))/length(ind)],...
    [4 4],[mean(VarAPV(ind))-std(VarAPV(ind))/length(ind) mean(VarAPV(ind))+std(VarAPV(ind))/length(ind)],...
    [7 7],[mean(VarPost(ind))-std(VarPost(ind))/length(ind) mean(VarPost(ind))+std(VarPost(ind))/length(ind)])

%%%%%%%%%%%
%%%%%%%%%



figure;hold on;
for i=1:12
subplot(6,2,i);hold on;
i=11
figure;hold on;
coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
plot(runningaverage(Experiment(i).timeACpre,20),runningaverage(coef*(Experiment(i).pitchACpre(round(median(Experiment(i).TargetingWN)),:)-mean(Experiment(i).pitchACpre(round(median(Experiment(i).TargetingWN)),:))),20),'*')
plot(runningaverage(Experiment(i).timeAPV,10),runningaverage(coef*(Experiment(i).pitchAPV(round(median(Experiment(i).TargetingWN)),:)-mean(Experiment(i).pitchACpre(round(median(Experiment(i).TargetingWN)),:))),10),'*','Color','g')
plot(runningaverage(Experiment(i).timeAPVwn,10),runningaverage(coef*(Experiment(i).pitchAPVwn(round(median(Experiment(i).TargetingWN)),:)-mean(Experiment(i).pitchACpre(round(median(Experiment(i).TargetingWN)),:))),10),'*','Color','r')
plot(runningaverage(Experiment(i).timeACpost,20),runningaverage(coef*(Experiment(i).pitchACpost(round(median(Experiment(i).TargetingWN)),:)-mean(Experiment(i).pitchACpre(round(median(Experiment(i).TargetingWN)),:))),20),'*','Color','k')
end
figure;hold on;
            subplot(312);hold on;
            base=[];
            ind=[1 2 4 5 6 7 8 10 12:24];
            for i=ind
                window=round(median(Experiment(i).TargetingWN));
                EtimeAPV=Experiment(i).timeAPVwn;
                [Esorted,sortedAPV]=sort(EtimeAPV);
                mpitchpre=mean(mean(Experiment(i).pitchACpre(window,:)));
                coef=1-2*(isequal(Experiment(i).DIR,'down')); % -1 if down, 1 if up
                a=Experiment(i).pitchAPV(window,round(end/2):end);
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                c=mpitchpre;
                base(i).data=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(1:round(end/3)))))
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                first(i)=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                b=mean(mean(Experiment(i).pitchACpre(window,:)));
                preAPV(i)=coef*(a-b);
                preAPVnorm(i)=preAPV(i)/c;
                a=mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(end/3):round(2*end/3)))));
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                second(i)=coef*(a-b)/c;
                %third(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(2*end/3):round(3*end/3)))))...
                   % -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                %fourth(i)=(coef*(mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(round(3*end/4):end))))...
                %    -mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)))))/mpitchpre;
                a=mean(mean(Experiment(i).pitchAPVwn(window,sortedAPV(end-15:end))));
                b=mean(mean(Experiment(i).pitchAPV(window,round(end/2):end)));
                finalpoint(i)=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchACpost(window,1:20)));
                b=mean(mean(Experiment(i).pitchACpre(window,:)));
                post1(i)=coef*(a-b)/c;
                a=mean(mean(Experiment(i).pitchACpost(window,1:postday(i,2)-1)));
                b=mean(mean(Experiment(i).pitchACpre(window,:)));
                post2(i)=coef*(a-b)/c;
                LearningNorm(ind)=Learning(ind)/c;
            end
            
     % 95th
     CIbase=[];
     for i=ind
            CIbase(i,:)=resampledCI(base(i).data,0.05,1);
     end
     CIbase=mean(CIbase(ind,:));
            CIfirst=resampledCI(first(ind),0.05,1);
            CIsecond=resampledCI(second(ind),0.05,1);
            %CIthird=resampledCI(third(ind),0.05);
            CIfinal=resampledCI(finalpoint(ind),0.05,1);
            CIpost1=resampledCI(post1(ind),0.05,1);
            CIpost2=resampledCI(post2(ind),0.05,1);
            
            
            figure;hold on;
            subplot(2,1,1);hold on;
            plot([1 1],[CIfirst],'r')
            plot([2 2],[CIsecond],'r')
            %plot([3 3],[CIthird],'r')
            plot([3 3],[CIfinal],'r')
            plot([4 4],[CIpost1],'r')
            plot([5 5],[CIpost2],'r')
            plot([0 1 2 3 4 5],[0 mean(first(ind)) mean(second(ind)) mean(finalpoint(ind)) mean(post1(ind)) mean(post2(ind))],'r')
            plot([-1 6],[0 0],'k')

        % 99th
%             CIfirst=resampledCI(first(ind),0.01,1);
%             CIsecond=resampledCI(second(ind),0.01,1);
%             %CIthird=resampledCI(third(ind),0.05);
%             CIfinal=resampledCI(finalpoint(ind),0.01,1);
%             CIpost1=resampledCI(post1(ind),0.01,1);
%             CIpost2=resampledCI(post2(ind),0.01,1);
            
            
            figure;hold on;
            plot([0 0],[CIbase],'r')
%             plot([1 1],[CIfirst],'r')
%             plot([2 2],[CIsecond],'r')
            %plot([3 3],[CIthird],'r')
            plot([3 3],[CIfinal],'r')
            plot([4 4],[CIpost1],'r')
            plot([5 5],[CIpost2],'r')
            plot([0 3 4 5],[0 mean(finalpoint(ind)) mean(post1(ind)) mean(post2(ind))],'r')
            plot([-1 6],[0 0],'k')

   basePC=[];
            
            for i=1:length(ExperimentPC)
                
                coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
                mpitchpre=mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)));
                basePC(i).data=(coef*(((ExperimentPC(i).pitchPre(round(median(ExperimentPC(i).time)),:)))-mpitchpre))/mpitchpre; 
                onethirdPC(i)=(coef*(mean(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,1:round(end/3))))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                twothirdsPC(i)=(coef*(mean(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,round(end/3):round(2*end/3))))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                finalpointPC(i)=(coef*(mean(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,end-15:end)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                postPC1(i)=(coef*(mean(mean(ExperimentPC(i).pitchPost(ExperimentPC(i).time,1:20)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;
                postPC2(i)=(coef*(mean(mean(ExperimentPC(i).pitchPost(ExperimentPC(i).time,:)))...
                    -mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)))))/mpitchpre;

            end
            for i=1:14
            pcCIbase(i,:)=resampledCI(basePC(i).data,0.1,1);
            end
            pcCIbase=mean(pcCIbase);
            pcCIfirst=resampledCI(onethirdPC,0.05,1);
            pcCIsecond=resampledCI(twothirdsPC,0.05,1);
            %CIthird=resampledCI(third(ind),0.05);
            pcCIfinal=resampledCI(finalpointPC,0.05,1);
            pcCIpost1=resampledCI(postPC1,0.05,1);
            pcCIpost2=resampledCI(postPC2,0.05,1);
            
            figure;hold on
%             plot([0.9 0.9],[pcCIfirst],'b')
%             plot([1.9 1.9],[pcCIsecond],'b')
            %plot([3 3],[CIthird],'r')
            plot([-0.1 -0.1],[pcCIbase],'b')
            plot([2.9 2.9],[pcCIfinal],'b')
            plot([3.9 3.9],[pcCIpost1],'b')
            plot([4.9 4.9],[pcCIpost2],'b')
            plot([-0.1 2.9 3.9 4.9],[0 mean(finalpointPC) mean(postPC1) mean(postPC2)],'b')

%%% one problem is that preAPV(ind) --- i.e. the offset effect --- is
%%% fairly large and in the adaptive direction and correlated with learning

% BUT - CONTROL NOTES



% 811 - figs
figure;hold on;
            plot([-0.1 -0.1],[pcCIbase],'b')
            plot([2.9 2.9],[pcCIfinal],'b')
            plot([3.4 3.4],[pcCIpost1],'b')
            plot([4.9 4.9],[pcCIpost2],'b')
            plot([-0.1 2.9 3.4 4.9],[0 mean(finalpointPC) mean(postPC1) mean(postPC2)],'b')
            plot([0 0],[CIbase],'r')
%             plot([1 1],[CIfirst],'r')
%             plot([2 2],[CIsecond],'r')
            %plot([3 3],[CIthird],'r')
            plot([3 3],[CIfinal],'r')
            plot([3.3 3.3],[CIpost1],'r')
            plot([5 5],[CIpost2],'r')
            plot([0 3 3.3 5],[0 mean(finalpoint(ind)) mean(post1(ind)) mean(post2(ind([1:6 8:end])))],'r')
            plot([-1 6],[0 0],'k')
figure;hold on;
bar([mean(VarPre(ind)) mean(VarAPV(ind)) mean(VarPost(ind))])
plot([1 1],[mean(VarPre(ind))-std(VarPre(ind))/length(ind) mean(VarPre(ind))+std(VarPre(ind))/length(ind)],...
    [2 2],[mean(VarAPV(ind))-std(VarAPV(ind))/length(ind) mean(VarAPV(ind))+std(VarAPV(ind))/length(ind)],...
    [3 3],[mean(VarPost(ind))-std(VarPost(ind))/length(ind) mean(VarPost(ind))+std(VarPost(ind))/length(ind)])

%%%%%% alternative


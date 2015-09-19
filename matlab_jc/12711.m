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
                
                
                
 % Figure 4A - Total hours (includes sleep)
for i=1:length(ExperimentPC)
    timepost=ExperimentPC(1).timePost-max(ExperimentPC(1).timeWN);
    postday
end
        subplot(221);figure;hold on
        for i=1:length(ExperimentPC)
%             colorchoice(i,:)=[rand rand rand];
            coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
            window=[ExperimentPC(i).time-30:ExperimentPC(i).time];
            % All data from the beginning of the next day onwards
            Etime=ExperimentPC(i).timePost;
            [x,sortedPost]=sort(Etime);
            Epitch=ExperimentPC(i).pitchPost(:,1:end);
            plot(runningaverage(Etime(sortedPost)-Etime(sortedPost(1)),10),...
                coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(ExperimentPC(i).pitchPre(window,:))),10)),'-','Color',colorchoice(i,:))
            timevals(i).data=Etime(sortedPost)-Etime(sortedPost(1));
            pitchvals(i).data=coef*(mean(Epitch(window,sortedPost))-mean(mean(ExperimentPC(i).pitchPre(window,:))));
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
            coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
            window=[ExperimentPC(i).on:ExperimentPC(i).off];
            Etime=ExperimentPC(i).timeACpost(:,postday(i,1):end);
            [Esorted,sortedPost]=sort(Etime);
          % Finds indices of gaps longer than 8hrs
            nights=find(diff(Esorted)>8);
          % Subtracts gaps
            for j=1:length(nights)
                gaplength=Esorted(nights(j)+1)-Esorted(nights(j));
                Esorted(nights(j)+1:end)=Esorted(nights(j)+1:end)-gaplength;
            end
          % plots points with appropriate x-axis
            Epitch=ExperimentPC(i).pitchACpost(:,postday(i,1):end);
            plot(runningaverage(Esorted-Esorted(1),20),...
                coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(ExperimentPC(i).pitchACpre(window,:))),20)),'Color',colorchoice(i,:))
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
                    coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
                    window=[ExperimentPC(i).on:ExperimentPC(i).off];
                    % All data from the beginning of the next day onwards
                    Etime=ExperimentPC(i).timeACpost(:,1:end);
                    [x,sortedPost]=sort(Etime);
                    Epitch=ExperimentPC(i).pitchACpost(:,1:end);
                    plot(runningaverage(Etime(sortedPost)-Etime(sortedPost(1)),20),...
                        coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(ExperimentPC(i).pitchACpre(window,:))),20)),'-','Color',colorchoice(i,:))
                    timevals(i).data=Etime(sortedPost)-Etime(sortedPost(1));
                    pitchvals(i).data=coef*(mean(Epitch(window,sortedPost))-mean(mean(ExperimentPC(i).pitchACpre(window,:))));
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
                    coef=1-2*(isequal(ExperimentPC(i).DIR,'down')); % -1 if down, 1 if up
                    window=[ExperimentPC(i).on:ExperimentPC(i).off];
                    Etime=ExperimentPC(i).timeACpost(:,1:end);
                    [Esorted,sortedPost]=sort(Etime);
                  % Finds indices of gaps longer than 8hrs
                    nights=find(diff(Esorted)>8);
                  % Subtracts gaps
                    for j=1:length(nights)
                        gaplength=Esorted(nights(j)+1)-Esorted(nights(j));
                        Esorted(nights(j)+1:end)=Esorted(nights(j)+1:end)-gaplength;
                    end
                  % plots points with appropriate x-axis
                    Epitch=ExperimentPC(i).pitchACpost(:,1:end);
                    plot(runningaverage(Esorted-Esorted(1),20),...
                        coef*(runningaverage(mean(Epitch(window,sortedPost))-mean(mean(ExperimentPC(i).pitchACpre(window,:))),20)),'Color',colorchoice(i,:))
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
% Divide into days

load /bulbul3/AnalysisMMAN/MANsummary3.mat

figure;hold on;
for m=1:length(Alldata)
    subplot(2,5,m);hold on;
    plot(Alldata(m).prelesionDays,Alldata(m).prelesionData,'k')
    ind1=find(Alldata(m).prelesionDays<0);
    ind2=find(Alldata(m).prelesionDays>=0 & Alldata(m).prelesionDays<=Alldata(m).postlesionWNendday);
    ind3=find(Alldata(m).prelesionDays>Alldata(m).prelesionWNendday);
    plot(Alldata(m).prelesionDays(ind1),Alldata(m).prelesionData(ind1),'b.','Markersize',20)
    plot(Alldata(m).prelesionDays(ind2),Alldata(m).prelesionData(ind2),'r.','Markersize',20)
    plot(Alldata(m).prelesionDays(ind3),Alldata(m).prelesionData(ind3),'b.','Markersize',20)
    if m~=4;xlim(Alldata(m).xlimpre);ylim(Alldata(m).ylimpre);end
    subplot(2,5,m+5);hold on;
    plot(Alldata(m).postlesionDays,Alldata(m).postlesionData,'k')    
    ind1=find(Alldata(m).postlesionDays<0);
    ind2=find(Alldata(m).postlesionDays>=0 & Alldata(m).postlesionDays<=Alldata(m).postlesionWNendday);
    ind3=find(Alldata(m).postlesionDays>Alldata(m).postlesionWNendday);
    plot(Alldata(m).postlesionDays(ind1),Alldata(m).postlesionData(ind1),'b.','Markersize',20)
    plot(Alldata(m).postlesionDays(ind2),Alldata(m).postlesionData(ind2),'r.','Markersize',20)
    plot(Alldata(m).postlesionDays(ind3),Alldata(m).postlesionData(ind3),'b.','Markersize',20)
    xlim(Alldata(m).xlimpost);ylim(Alldata(m).ylimpost)
end


figure;hold on;
% TW data
        subplot(2,5,1);hold on;
        plot([1:4],sqout(11).BASALL,'.','Markersize',20)
        plot(5:8,sqout(11).WNVALS,'r.','Markersize',20)
        plot(10,sqout(11).REC,'.','Markersize',20)
        plot([1:8 10],[sqout(11).BASALL; sqout(11).WNVALS; sqout(11).REC],'k')

        subplot(2,5,2);hold on;
        plot([1:2],sqout(14).BASALL,'.','Markersize',20)
        plot(7:8,sqout(14).WNVALS,'r.','Markersize',20)
        plot(11,sqout(14).REC,'.','Markersize',20)
        plot([1:2 7:8 11],[sqout(14).BASALL; sqout(14).WNVALS; sqout(14).REC],'k')

% Branch points
        % Pre-lesion
            for i=[4]
                subplot(2,5,i);hold on;
                % Branch points
                clear days relday mnprct
                    days{3}=[1 find(diff(Exp(i).SYNpre{3,1})>8) length(Exp(i).SYNpre{3,1})];
                    days{2}=[1 find(diff(Exp(i).SYNpre{2,1})>8) length(Exp(i).SYNpre{2,1})];
                    days{1}=[1 find(diff(Exp(i).SYNpre{1,1})>8) length(Exp(i).SYNpre{1,1})];
                    for j=1:3
                        count=0;
                        for k=1:length(days{j})-1
                            if days{j}(k+1)-days{j}(k)>20
                                count=count+1;
                                relday{j}(count)=floor((mean(Exp(i).SYNpre{j,1}(days{j}(k):days{j}(k+1)))-min(Exp(i).SYNpre{2,1}))/24);
                                mnprct{j}(count)=mean(Exp(i).SYNpre{j,2}(days{j}(k):days{j}(k+1)));
                            end
                        end
                    end
                    plot(relday{1},mnprct{1},'.','Markersize',20)
                    plot(relday{2},mnprct{2},'r.','Markersize',20)
                    plot(relday{3},mnprct{3},'.','Markersize',20)    
                    plot([relday{1} relday{2} relday{3}],[mnprct{1} mnprct{2} mnprct{3}],'k')
              
            end


        % Post-lesion
            for i=[1 3 4]
                subplot(2,5,5+i);hold on;
                % Branch points
                clear days relday mnprct
                    days{3}=[1 find(diff(Exp(i).SYNpost{3,1})>8) length(Exp(i).SYNpost{3,1})];
                    days{2}=[1 find(diff(Exp(i).SYNpost{2,1})>8) length(Exp(i).SYNpost{2,1})];
                    days{1}=[1 find(diff(Exp(i).SYNpost{1,1})>8) length(Exp(i).SYNpost{1,1})];
                    for j=1:3
                        count=0;
                        for k=1:length(days{j})-1
                            if days{j}(k+1)-days{j}(k)>20
                                count=count+1;
                                relday{j}(count)=floor((mean(Exp(i).SYNpost{j,1}(days{j}(k):days{j}(k+1)))-min(Exp(i).SYNpost{2,1}))/24);
                                mnprct{j}(count)=mean(Exp(i).SYNpost{j,2}(days{j}(k):days{j}(k+1)));
                            end
                        end
                    end
                    plot(relday{1},mnprct{1},'.','Markersize',20)
                    plot(relday{2},mnprct{2},'r.','Markersize',20)
                    plot(relday{3},mnprct{3},'.','Markersize',20)       
                    plot([relday{1} relday{2} relday{3}],[mnprct{1} mnprct{2} mnprct{3}],'k')
           
            end

    % Repeats
% Pre-lesion
        for i=[1]
                subplot(2,5,5);hold on;
                clear days relday mnrep
                    days{3}=[1 find(diff(Exp(i).REPpre{3,1})>8) length(Exp(i).REPpre{3,1})];
                    days{2}=[1 find(diff(Exp(i).REPpre{2,1})>8) length(Exp(i).REPpre{2,1})];
                    days{1}=[1 find(diff(Exp(i).REPpre{1,1})>8) length(Exp(i).REPpre{1,1})];
                    for j=1:3
                        count=0;
                        for k=1:length(days{j})-1
                            if days{j}(k+1)-days{j}(k)>20
                                count=count+1;
                                relday{j}(count)=floor((mean(Exp(i).REPpre{j,1}(days{j}(k):days{j}(k+1)))-min(Exp(i).REPpre{2,1}))/24);
                                mnrep{j}(count)=mean(Exp(i).REPpre{j,2}(days{j}(k):days{j}(k+1)));
                            end
                        end
                    end
                    plot(relday{1},mnrep{1},'.','Markersize',20)
                    plot(relday{2},mnrep{2},'r.','Markersize',20)
                    plot(relday{3},mnrep{3},'.','Markersize',20)   
                    plot([relday{1} relday{2} relday{3}],[mnrep{1} mnrep{2} mnrep{3}],'k')
           
            end


        % Post-lesion
            for i=[1 2]
                if i==2
                    subplot(2,5,5+i);hold on;
                end
                if i==1
                    subplot(2,5,10);hold on;
                end
                % Branch points
                clear days relday mnrep
                    days{3}=[1 find(diff(Exp(i).REPpost{3,1})>8) length(Exp(i).REPpost{3,1})];
                    days{2}=[1 find(diff(Exp(i).REPpost{2,1})>8) length(Exp(i).REPpost{2,1})];
                    days{1}=[1 find(diff(Exp(i).REPpost{1,1})>8) length(Exp(i).REPpost{1,1})];
                    for j=1:3
                        count=0;
                        for k=1:length(days{j})-1
                            if days{j}(k+1)-days{j}(k)>20
                                count=count+1;
                                relday{j}(count)=floor((mean(Exp(i).REPpost{j,1}(days{j}(k):days{j}(k+1)))-min(Exp(i).REPpost{2,1}))/24);
                                mnrep{j}(count)=mean(Exp(i).REPpost{j,2}(days{j}(k):days{j}(k+1)));
                            end
                        end
                    end
                    plot(relday{1},mnrep{1},'.','Markersize',20)
                    plot(relday{2},mnrep{2},'r.','Markersize',20)
                    plot(relday{3},mnrep{3},'.','Markersize',20)  
                    plot([relday{1} relday{2} relday{3}],[mnrep{1} mnrep{2} mnrep{3}],'k')
            end
    
    % FF


    
    
    
    
    
    
    
    
%%%%%% bk2bk10
i=4;
figure; hold on;
% 
subplot(2,3,2);hold on;
    ravgsyn=60;
    plot(runningaverage(Exp(i).SYNpre{1,1}-min(Exp(i).SYNpre{2,1}),ravgsyn),runningaverage(Exp(i).SYNpre{1,2},ravgsyn),'b')
    plot(runningaverage(Exp(i).SYNpre{2,1}-min(Exp(i).SYNpre{2,1}),ravgsyn),runningaverage(Exp(i).SYNpre{2,2},ravgsyn),'r')
    plot(runningaverage(Exp(i).SYNpre{3,1}-min(Exp(i).SYNpre{2,1}),ravgsyn),runningaverage(Exp(i).SYNpre{3,2},ravgsyn),'k')


subplot(2,3,5);hold on;
    ravgsyn=60;
    plot(runningaverage(Exp(i).SYNpost{1,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{1,2},ravgsyn),'b')
    plot(runningaverage(Exp(i).SYNpost{2,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{2,2},ravgsyn),'r')
    plot(runningaverage(Exp(i).SYNpost{3,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{3,2},ravgsyn),'k')


%%%%%% b65o7
i=1;
figure;hold on;
% FF
subplot(2,3,4);hold on;
    ravgFF=40;
    plot(runningaverage(timing3(Exp(i).FFpost{1,1})-min(timing3(Exp(i).FFpost{2,1})),ravgFF),...
        runningmedian(mean(Exp(i).FFpost{1,2}(Exp(i).FFwindow,:)),ravgFF),'b')
    plot(runningaverage(timing3(Exp(i).FFpost{2,1})-min(timing3(Exp(i).FFpost{2,1})),ravgFF),...
        runningmedian(mean(Exp(i).FFpost{2,2}(Exp(i).FFwindow,:)),ravgFF),'r')
    plot(runningaverage(timing3(Exp(i).FFpost{3,1})-min(timing3(Exp(i).FFpost{2,1})),ravgFF),...
        runningmedian(mean(Exp(i).FFpost{3,2}(Exp(i).FFwindow,:)),ravgFF),'k')

% Syntax
subplot(2,3,5);hold on;
    ravgsyn=40;
    plot(runningaverage(Exp(i).SYNpost{1,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{1,2},ravgsyn),'b')
    plot(runningaverage(Exp(i).SYNpost{2,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{2,2},ravgsyn),'r')
    plot(runningaverage(Exp(i).SYNpost{3,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{3,2},ravgsyn),'k')

% Repeats
subplot(2,3,3);hold on;
    ravgrep=40;
    plot(runningaverage(Exp(i).REPpre{1,1}-min(Exp(i).REPpre{2,1}),ravgrep),runningaverage(Exp(i).REPpre{1,2},ravgrep),'b')
    plot(runningaverage(Exp(i).REPpre{2,1}-min(Exp(i).REPpre{2,1}),ravgrep),runningaverage(Exp(i).REPpre{2,2},ravgrep),'r')
    plot(runningaverage(Exp(i).REPpre{3,1}-min(Exp(i).REPpre{2,1}),ravgrep),runningaverage(Exp(i).REPpre{3,2},ravgrep),'k')
subplot(2,3,6);hold on;
    ravgrep=40;
    plot(runningaverage(Exp(i).REPpost{1,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpost{1,2},ravgrep),'b')
    plot(runningaverage(Exp(i).REPpost{2,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpost{2,2},ravgrep),'r')
    plot(runningaverage(Exp(i).REPpost{3,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpost{3,2},ravgrep),'k')

%%%%%% w87o32
i=2; 
figure;hold on;
ravgrep=40;
    plot(runningaverage(Exp(i).REPpre{1,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpre{1,2},ravgrep),'b')
    plot(runningaverage(Exp(i).REPpost{1,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpost{1,2},ravgrep),'b')
    plot(runningaverage(Exp(i).REPpost{2,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpost{2,2},ravgrep),'r')
    plot(runningaverage(Exp(i).REPpost{3,1}-min(Exp(i).REPpost{2,1}),ravgrep),runningaverage(Exp(i).REPpost{3,2},ravgrep),'k')


%%%%%% pu19bk81
i=3
figure;hold on;
% FF
subplot(2,3,4);hold on;
    ravgFF=40;
    plot(runningaverage(timing3(Exp(i).FFpost{1,1})-min(timing3(Exp(i).FFpost{2,1})),ravgFF),...
        runningmedian(mean(Exp(i).FFpost{1,2}(Exp(i).FFwindow,:)),ravgFF),'b')
    plot(runningaverage(timing3(Exp(i).FFpost{2,1})-min(timing3(Exp(i).FFpost{2,1})),ravgFF),...
        runningmedian(mean(Exp(i).FFpost{2,2}(Exp(i).FFwindow,:)),ravgFF),'r')
    plot(runningaverage(timing3(Exp(i).FFpost{3,1})-min(timing3(Exp(i).FFpost{2,1})),ravgFF),...
        runningmedian(mean(Exp(i).FFpost{3,2}(Exp(i).FFwindow,:)),ravgFF),'k')

% Syntax
subplot(2,3,5);hold on;
    ravgsyn=40;
    plot(runningaverage(Exp(i).SYNpost{1,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{1,2},ravgsyn),'b')
    plot(runningaverage(Exp(i).SYNpost{2,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{2,2},ravgsyn),'r')
    plot(runningaverage(Exp(i).SYNpost{3,1}-min(Exp(i).SYNpost{2,1}),ravgsyn),runningaverage(Exp(i).SYNpost{3,2},ravgsyn),'k')




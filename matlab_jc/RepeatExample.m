
    %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% repeat length example
  %%% July 16, 2011
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  
       clear all
    load /bulbul3/AnalysisRepeats/RepeatExampleReplength.mat
           

mnpost1=mean(distntPOST(1:105));
sempost1=std(distntPOST(1:105))/sqrt(105);
mnpost2=mean(distntPOST(106:267));
sempost2=std(distntPOST(106:267))/sqrt(162);
mnpost3=mean(distntPOST(268:404));
sempost3=std(distntPOST(268:404))/sqrt(137);
mnpost4=mean(distntPOST2(7:end));
sempost4=std(distntPOST2(7:end))/sqrt(241);

                    figure;hold on;
                    i=7;
                        alldays=[];
                        plot(Exp1(i).predays(5:7),Exp1(i).replong(Exp1(i).predays(5:7)),'.','Markersize',25)
                        plot([Exp1(i).predays(5:7);Exp1(i).predays(5:7)],...
                            [Exp1(i).replong(Exp1(i).predays(5:7))-Exp1(i).semreplong(Exp1(i).predays(5:7));Exp1(i).replong(Exp1(i).predays(5:7))+Exp1(i).semreplong(Exp1(i).predays(5:7))],'k-')    
                        alldays=[alldays Exp1(i).predays(5:7) Exp1(i).wndays 15.1 15.8];
                        plot(Exp1(i).wndays,Exp1(i).replong(Exp1(i).wndays),'r.','Markersize',25)
                        plot([Exp1(i).wndays;Exp1(i).wndays],...
                            [Exp1(i).replong(Exp1(i).wndays)-Exp1(i).semreplong(Exp1(i).wndays);Exp1(i).replong(Exp1(i).wndays)+Exp1(i).semreplong(Exp1(i).wndays)],'k-')    
                        plot([15.1 15.8 ],[mnpost1 mnpost2 ],'.','Markersize',25)
                        plot([15.1  15.8 ;15.1 15.8 ],...
                            [mnpost1-sempost1 mnpost2-sempost2 ;mnpost1+sempost1 mnpost2+sempost2],'k-')    
                        plot(alldays,[Exp1(i).replong(Exp1(i).predays(5:7)) Exp1(i).replong(Exp1(i).wndays) mnpost1 mnpost2 ],'k','Linewidth',2)
                        ylim([2 5])
                        xlim([7 18])
                        plot([7 18],[mean(Exp1(i).replong(Exp1(i).predays(5:7))) mean(Exp1(i).replong(Exp1(i).predays(5:7)))],':')

                    % HISTOGRAMS
                        clear d
                        daycount=0;
                        for j=1:3;
                        if Exp(i).JC==0 % Time data for Evren vs. JC
                            unidays=find(diff(Exp(i).time{j})>0.3); % Evren time
                            nightlength=1/3;
                        else
                            unidays=find(diff(Exp(i).time{j})>8); % JC time
                            nightlength=8;
                        end
                        unistart=[1 unidays+1];
                        uniend=[unidays length(Exp(i).time{j})];
                        for ii=1:length(unistart) % For each day
                            daycount=1+ceil((min(Exp(i).time{j}(unistart(ii):uniend(ii)))-min(Exp(i).time{1}))/(3*nightlength));
                            allrpdata{daycount}=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                            group(daycount)=j;
                        end
                        end
                        % For each day, figure out all parameters
                        for j=1:length(allrpdata)
                            [b,a]=hist(allrpdata{j},[0:1:15]);
                            [c,d(j,:)]=stairs(a,b/sum(b));
                            mnrpvls(j)=mean(allrpdata{j});
                            sdrpvls(j)=std(allrpdata{j});
                        end
                        figure;hold on;
                        plot(c,d([8 10],:),'Linewidth',2,'Color','b')
                        plot(c,d(14,:),'Linewidth',4,'Color','r')
                        plot(c,d(16,:),'Linewidth',2,'Color','k')    
                        %plot(c,d(find(group==2),:),'Linewidth',2,'Color','r')
                        %plot(c,d(find(group==3),:),'Linewidth',2,'Color','k')
                        xlim([0 12])
                        ylim([0 0.9])

                        figure;hold on;
                        plot(c,d([1:3 8:10],:),'Linewidth',2)
                        xlim([0 12])
                        ylim([0 0.45])
                        plot(mnrpvls([1:3 8:10]),0.35,'.','Markersize',25)

                    % TRANSITION PROBABILITIES
                  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% ptrans example
  %%% July 16, 2011
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  clear all
  load /bulbul3/AnalysisRepeats/RepeatExamplePtrans.mat
  
                % Do mean pre and mean post to simplify    
                    i=10;
                    preMN=mean(Exp1(i).ptrans([2 4],1:6));
                  % PRE
             
                       rpdata=Exp(10).rplength{1};
                       day=1;
                       for k=1:1:20 % for each transition
                            ptrans0(day,k)=length(find(rpdata>k))/length(find(rpdata>k-1));
                            samplesizeptrans(day,k)=length(find(rpdata>k));
                            clear ptransBS
                            for bs=1:1000
                                scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                                scrambleRPL=rpdata(scrambleIND);
                                ptransBS(bs)=length(find(scrambleRPL>k))/length(find(scrambleRPL>k-1));
                            end
                            semptrans0(day,k)=std(ptransBS);
                       end

               % WN day 2 afternoon --- Exp(10).rplength{2}(84:111);
                       rpdata=Exp(10).rplength{2}(84:111);
                       day=1;
                       for k=1:1:20 % for each transition
                            ptrans1(day,k)=length(find(rpdata>k))/length(find(rpdata>k-1));
                            samplesizeptrans(day,k)=length(find(rpdata>k));
                            clear ptransBS
                            for bs=1:1000
                                scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                                scrambleRPL=rpdata(scrambleIND);
                                ptransBS(bs)=length(find(scrambleRPL>k))/length(find(scrambleRPL>k-1));
                            end
                            semptrans1(day,k)=std(ptransBS);
                       end
               % Post day first hours --- Exp(10).rplength{2}(84:111);
                       rpdata=distntPOST(1:52);
                       day=1;
                       for k=1:1:20 % for each transition
                            ptrans2(day,k)=length(find(rpdata>k))/length(find(rpdata>k-1));
                            samplesizeptrans(day,k)=length(find(rpdata>k));
                            clear ptransBS
                            for bs=1:1000
                                scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                                scrambleRPL=rpdata(scrambleIND);
                                ptransBS(bs)=length(find(scrambleRPL>k))/length(find(scrambleRPL>k-1));
                            end
                            semptrans2(day,k)=std(ptransBS);
                       end
                    
                    wnday1MN=Exp1(i).ptrans(5,1:6);
                    postMN=(Exp1(i).ptrans([8],1:6));
                    
                    allptrans=[ptrans0(1:6);wnday1MN(1:6);ptrans1(1:6);ptrans2(1:6);postMN(1:6)];
                    allsemptrans=[semptrans0(1:6);Exp1(i).semptrans(5,1:6);semptrans1(1:6);semptrans2(1:6);Exp1(i).semptrans(8,1:6)];
                    
                    figure;hold on;
                    plot([1 2 3 4 5],allptrans,'Linewidth',2)
                    k=1;plot([[1 2 3 4 5];[1 2 3 4 5]],[allptrans([1 2 3 4 5],k)-allsemptrans([1 2 3 4 5],k),allptrans([1 2 3 4 5],k)+allsemptrans([1 2 3 4 5],k)]','b')
                    k=2;plot([[1 2 3 4 5];[1 2 3 4 5]],[allptrans([1 2 3 4 5],k)-allsemptrans([1 2 3 4 5],k),allptrans([1 2 3 4 5],k)+allsemptrans([1 2 3 4 5],k)]','g')
                    k=3;plot([[1 2 3 4 5];[1 2 3 4 5]],[allptrans([1 2 3 4 5],k)-allsemptrans([1 2 3 4 5],k),allptrans([1 2 3 4 5],k)+allsemptrans([1 2 3 4 5],k)]','r')
                    k=4;plot([[1 2 3 4 5];[1 2 3 4 5]],[allptrans([1 2 3 4 5],k)-allsemptrans([1 2 3 4 5],k),allptrans([1 2 3 4 5],k)+allsemptrans([1 2 3 4 5],k)]','c')
                    k=5;plot([[1 2 3 4 5];[1 2 3 4 5]],[allptrans([1 2 3 4 5],k)-allsemptrans([1 2 3 4 5],k),allptrans([1 2 3 4 5],k)+allsemptrans([1 2 3 4 5],k)]','m')
                    k=6;plot([[1 2 3 4 5];[1 2 3 4 5]],[allptrans([1 2 3 4 5],k)-allsemptrans([1 2 3 4 5],k),allptrans([1 2 3 4 5],k)+allsemptrans([1 2 3 4 5],k)]','y')
                    plot([1 2 3 4 5],allptrans([1 2 3 4 5],1:6),'.','Markersize',25)
                    plot([1;5],[(allptrans(1,1:6));(allptrans(1,1:6))],':')
                        ylim([0 1.1])
                        xlim([0 6])
                        figure;plot(allptrans(1,:))
                %%%% Histo
                                            clear d
                        daycount=0;
                        for j=1:3;
                        if Exp(i).JC==0 % Time data for Evren vs. JC
                            unidays=find(diff(Exp(i).time{j})>0.3); % Evren time
                            nightlength=1/3;
                        else
                            unidays=find(diff(Exp(i).time{j})>8); % JC time
                            nightlength=8;
                        end
                        unistart=[1 unidays+1];
                        uniend=[unidays length(Exp(i).time{j})];
                        for ii=1:length(unistart) % For each day
                            daycount=1+ceil((min(Exp(i).time{j}(unistart(ii):uniend(ii)))-min(Exp(i).time{1}))/(3*nightlength));
                            allrpdata{daycount}=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                            group(daycount)=j;
                        end
                        end
                        % For each day, figure out all parameters
                        for j=1:length(allrpdata)
                            [b,a]=hist(allrpdata{j},[0:1:15]);
                            [c,d(j,:)]=stairs(a,b/sum(b));
                        end
                        figure;hold on;
                        plot(c,d([1 ],:),'Linewidth',2,'Color','b')
%                         plot(c,d(6,:),'Linewidth',4,'Color','r')
%                         plot(c,d(9,:),'Linewidth',2,'Color','k')    
                        %plot(c,d(find(group==2),:),'Linewidth',2,'Color','r')
                        %plot(c,d(find(group==3),:),'Linewidth',2,'Color','k')
                        xlim([0 8.9])
                        ylim([0 0.4])

                        figure;hold on;
                        plot(c,d([1 3 4],:),'Linewidth',2)
                        xlim([0 12])
                        ylim([0 0.6])

    
 % example 19 - b65o7

       figure;hold on;
                    
%                     i=7;
%                         alldays=[];
%                         plot(Exp1(i).predays(2:7),Exp1(i).replong(Exp1(i).predays(2:7)),'.','Markersize',25)
%                         plot([Exp1(i).predays(2:7);Exp1(i).predays(2:7)],...
%                             [Exp1(i).replong(Exp1(i).predays(2:7))-Exp1(i).semreplong(Exp1(i).predays(2:7));Exp1(i).replong(Exp1(i).predays(2:7))+Exp1(i).semreplong(Exp1(i).predays(2:7))],'k-')    
%                         plot(Exp1(i).predays(2:7),Exp1(i).replong(Exp1(i).predays(2:7)))
%                         ylim([2 5])
%                         xlim([1 11])
% 
% 

                    figure;hold on;
                    i=10;
                        alldays=[];
                        plot(Exp1(i).predays,Exp1(i).replong(Exp1(i).predays),'.','Markersize',25)
                        plot([Exp1(i).predays;Exp1(i).predays],...
                            [Exp1(i).replong(Exp1(i).predays)-Exp1(i).semreplong(Exp1(i).predays);Exp1(i).replong(Exp1(i).predays)+Exp1(i).semreplong(Exp1(i).predays)],'k-')    
                        alldays=[alldays Exp1(i).predays Exp1(i).wndays Exp1(i).postdays];
                        plot(Exp1(i).wndays,Exp1(i).replong(Exp1(i).wndays),'r.','Markersize',25)
                        plot([Exp1(i).wndays;Exp1(i).wndays],...
                            [Exp1(i).replong(Exp1(i).wndays)-Exp1(i).semreplong(Exp1(i).wndays);Exp1(i).replong(Exp1(i).wndays)+Exp1(i).semreplong(Exp1(i).wndays)],'k-')    
                        plot(Exp1(i).postdays,Exp1(i).replong(Exp1(i).postdays),'.','Markersize',25)
                        plot([Exp1(i).postdays;Exp1(i).postdays],...
                            [Exp1(i).replong(Exp1(i).postdays)-Exp1(i).semreplong(Exp1(i).postdays);Exp1(i).replong(Exp1(i).postdays)+Exp1(i).semreplong(Exp1(i).postdays)],'k-')    
                        plot(alldays,Exp1(i).replong(alldays),'k','Linewidth',2)
                    plot(runningaverage(Exp(i).time{3}/24,20)-min(Exp(i).time{3}/24)+6.5,runningaverage(Exp(i).rplength{3},20))
    
                        ylim([2 5])
                        xlim([7 18])
                        plot([7 18],[mean(Exp1(i).replong(Exp1(i).predays)) mean(Exp1(i).replong(Exp1(i).predays))],':')

                    % HISTOGRAMS
                        clear d
                        daycount=0;
                        for j=1:3;
                        if Exp(i).JC==0 % Time data for Evren vs. JC
                            unidays=find(diff(Exp(i).time{j})>0.3); % Evren time
                            nightlength=1/3;
                        else
                            unidays=find(diff(Exp(i).time{j})>8); % JC time
                            nightlength=8;
                        end
                        unistart=[1 unidays+1];
                        uniend=[unidays length(Exp(i).time{j})];
                        for ii=1:length(unistart) % For each day
                            daycount=1+ceil((min(Exp(i).time{j}(unistart(ii):uniend(ii)))-min(Exp(i).time{1}))/(3*nightlength));
                            allrpdata{daycount}=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                            group(daycount)=j;
                        end
                        end
                        % For each day, figure out all parameters
                        for j=1:length(allrpdata)
                            [b,a]=hist(allrpdata{j},[0:1:15]);
                            [c,d(j,:)]=stairs(a,b/sum(b));
                        end
                        figure;hold on;
                        plot(c,d([1 3 4],:),'Linewidth',2,'Color','b')
                        plot(c,d(6,:),'Linewidth',4,'Color','r')
                        plot(c,d(9,:),'Linewidth',2,'Color','k')    
                        %plot(c,d(find(group==2),:),'Linewidth',2,'Color','r')
                        %plot(c,d(find(group==3),:),'Linewidth',2,'Color','k')
                        xlim([0 12])
                        ylim([0 0.9])

                        figure;hold on;
                        plot(c,d([1 3 4],:),'Linewidth',2)
                        xlim([0 12])
                        ylim([0 0.6])

                    % TRANSITION PROBABILITIES
                      figure;hold on;                       
                        

%%%%%%
%%%%%%
% Example figure 
            load /bulbul1/LMAN_microstimulation/Data020912.mat
        ExpIndex=4;
            figure;hold on;
            catchfrac1=5;catchfrac2=1;ravgwin=50;
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
                        inds=[startpoints(j):catchfrac1:endpoints(j)];
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
                        inds=[startpoints(j):catchfrac2:endpoints(j)];           
                        plot(timedata(inds),FFdata(inds),'.','markersize',7,'Color','r')
                    end
                end
            end
            
         plot(runningaverage(timing3(IntExp(4).fvalsCTlearning),ravgwin),runningaverage(IntExp(4).FFdataCTlearning,ravgwin),'b','Linewidth',3)
         plot(runningaverage(timing3(IntExp(4).fvalsCTbaseline),ravgwin),runningaverage(IntExp(4).FFdataCTbaseline,ravgwin),'b','Linewidth',3)
         plot(runningaverage(timing3(IntExp(4).fvalsFBlearning),ravgwin),runningaverage(IntExp(4).FFdataFBlearning,ravgwin),'r','Linewidth',3)
         plot(runningaverage(timing3(IntExp(4).fvalsFBbaseline),ravgwin),runningaverage(IntExp(4).FFdataFBbaseline,ravgwin),'r','Linewidth',3)
         
         
            %%%%%%%%%%%%%%%%
            for i=1:length(IntExp(ExpIndex).FFmorningfitFBbaseline)
                t1=timing4(IntExp(ExpIndex).AllDaysFBbaseline(i).fv(1));
                t2=timing4(IntExp(ExpIndex).AllDaysFBbaseline(i).fv(end));
                plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBbaseline(i);IntExp(ExpIndex).FFeveningfitFBbaseline(i)],'r-','Linewidth',4)
            end

            for i=1:length(IntExp(ExpIndex).FFmorningfitFBlearning)
                t1=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(1));
                if i>1;plot([t2;t1],[IntExp(ExpIndex).FFeveningfitFBlearning(i-1);IntExp(ExpIndex).FFmorningfitFBlearning(i)],'k-','Linewidth',3);end             
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
                if i>1;plot([t2;t1],[IntExp(ExpIndex).FFeveningfitCTlearning(i-1);IntExp(ExpIndex).FFmorningfitCTlearning(i)],'k-','Linewidth',3);end
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
            for i=2%1:5
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
            
            learnday/(learnday+learnnight)
            transferday/(transferday+transfernight)
            
% Summary figure

%%%%%%%%%%%%%%
% "Dynamics of consolidation"
%%%%%%%%%%%%%%%
% Pretty straightforward, load the relevant file, daynight.m does the
% calculations very quickly, then plot the running avg'd data and linear
% fits for each day.
%%%%%%%%%%%%%%
% Things to keep in mind:
% 1. quality of reversion, since poor quality inactivations would be confounded by
% dynamics in the AFP component of learning (as opposed to dynamics of
% consolidation).
% 2. baseline dynamics
%%%%%%%%%%%%%%
clear all
load /bulbul1/IntExp030612.mat

%TW data - already processed
load /bulbul1/r87g80rvsd.mat


for i=1:length(vls)
vls(i,6)=fvpt_com(i).TRIG(1);
end
IntExp(5).bird='r87g80';IntExp(6).bird='r87g80';IntExp(7).bird='r87g80';IntExp(8).bird='r87g80';

IntExp(5).dir=1;IntExp(6).dir=-1;IntExp(7).dir=1;IntExp(8).dir=1;

CTind=crct_com(1:53246);
FBind=crfb_com(1:8919);
vlsCT=vls(CTind,:);
vlsFB=vls(FBind,:);
vlsCTexp1=vlsCT(find(vlsCT(:,1)<734500),:);
vlsFBexp1=vlsFB(find(vlsFB(:,1)<734500),:);
vlsCTexp2=vlsCT(find(vlsCT(:,1)>734500),:);
vlsFBexp2=vlsFB(find(vlsFB(:,1)>734500),:);
[IntExp(5).deltadayFBlearning,IntExp(5).deltanightFBlearning,IntExp(5).AllDaysFBlearning,IntExp(5).FFmorningfitFBlearning,IntExp(5).FFeveningfitFBlearning]...
    =daynightTW(vlsFBexp1);
[IntExp(6).deltadayFBlearning,IntExp(6).deltanightFBlearning,IntExp(6).AllDaysFBlearning,IntExp(6).FFmorningfitFBlearning,IntExp(6).FFeveningfitFBlearning]...
    =daynightTW(vlsFBexp2);
[IntExp(5).deltadayCTlearning,IntExp(5).deltanightCTlearning,IntExp(5).AllDaysCTlearning,IntExp(5).FFmorningfitCTlearning,IntExp(5).FFeveningfitCTlearning]...
    =daynightTW(vlsCTexp1);
[IntExp(6).deltadayCTlearning,IntExp(6).deltanightCTlearning,IntExp(6).AllDaysCTlearning,IntExp(6).FFmorningfitCTlearning,IntExp(6).FFeveningfitCTlearning]...
    =daynightTW(vlsCTexp2);

IntExp(7).deltadayFBlearning=IntExp(5).deltadayFBlearning(13:14);IntExp(7).deltadayFBbaseline=IntExp(5).deltadayFBlearning(11:2);
IntExp(7).deltanightFBlearning=IntExp(5).deltanightFBlearning(13:14);IntExp(7).deltanightFBbaseline=IntExp(5).deltanightFBlearning(11:12);
IntExp(7).FFmorningfitFBlearning=IntExp(5).FFmorningfitFBlearning(13:14);IntExp(7).FFmorningfitFBbaseline=IntExp(5).FFmorningfitFBlearning(11:12);
IntExp(7).AllDaysFBlearning=IntExp(5).AllDaysFBlearning(13:14);IntExp(7).AllDaysFBbaseline=IntExp(5).AllDaysFBlearning(11:12);
IntExp(7).FFeveningfitFBlearning=IntExp(5).FFeveningfitFBlearning(13:14);IntExp(7).FFeveningfitFBbaseline=IntExp(5).FFeveningfitFBlearning(11:12);
IntExp(7).deltadayCTlearning=IntExp(5).deltadayCTlearning(13:14);IntExp(7).deltadayCTbaseline=IntExp(5).deltadayCTlearning(11:12);
IntExp(7).deltanightCTlearning=IntExp(5).deltanightCTlearning(13:14);IntExp(7).deltanightCTbaseline=IntExp(5).deltanightCTlearning(11:12);
IntExp(7).FFmorningfitCTlearning=IntExp(5).FFmorningfitCTlearning(13:14);IntExp(7).FFmorningfitCTbaseline=IntExp(5).FFmorningfitCTlearning(11:12);
IntExp(7).AllDaysCTlearning=IntExp(5).AllDaysCTlearning(13:14);IntExp(7).AllDaysCTbaseline=IntExp(5).AllDaysCTlearning(11:12);
IntExp(7).FFeveningfitCTlearning=IntExp(5).FFeveningfitCTlearning(13:14);IntExp(7).FFeveningfitCTbaseline=IntExp(5).FFeveningfitCTlearning(11:12);

IntExp(5).deltadayFBlearning=IntExp(5).deltadayFBlearning(4:10);IntExp(5).deltadayFBbaseline=IntExp(5).deltadayFBlearning(1:2);
IntExp(5).deltanightFBlearning=IntExp(5).deltanightFBlearning(4:10);IntExp(5).deltanightFBbaseline=IntExp(5).deltanightFBlearning(1:2);
IntExp(5).FFmorningfitFBlearning=IntExp(5).FFmorningfitFBlearning(4:10);IntExp(5).FFmorningfitFBbaseline=IntExp(5).FFmorningfitFBlearning(1:2);
IntExp(5).AllDaysFBlearning=IntExp(5).AllDaysFBlearning(4:10);IntExp(5).AllDaysFBbaseline=IntExp(5).AllDaysFBlearning(1:2);
IntExp(5).FFeveningfitFBlearning=IntExp(5).FFeveningfitFBlearning(4:10);IntExp(5).FFeveningfitFBbaseline=IntExp(5).FFeveningfitFBlearning(1:2);
IntExp(5).deltadayCTlearning=IntExp(5).deltadayCTlearning(4:10);IntExp(5).deltadayCTbaseline=IntExp(5).deltadayCTlearning(1:2);
IntExp(5).deltanightCTlearning=IntExp(5).deltanightCTlearning(4:10);IntExp(5).deltanightCTbaseline=IntExp(5).deltanightCTlearning(1:2);
IntExp(5).FFmorningfitCTlearning=IntExp(5).FFmorningfitCTlearning(4:10);IntExp(5).FFmorningfitCTbaseline=IntExp(5).FFmorningfitCTlearning(1:2);
IntExp(5).AllDaysCTlearning=IntExp(5).AllDaysCTlearning(4:10);IntExp(5).AllDaysCTbaseline=IntExp(5).AllDaysCTlearning(1:2);
IntExp(5).FFeveningfitCTlearning=IntExp(5).FFeveningfitCTlearning(4:10);IntExp(5).FFeveningfitCTbaseline=IntExp(5).FFeveningfitCTlearning(1:2);

IntExp(8).deltadayFBlearning=IntExp(6).deltadayFBlearning(16:19);IntExp(8).deltadayFBbaseline=IntExp(6).deltadayFBlearning(14:15);
IntExp(8).deltanightFBlearning=IntExp(6).deltanightFBlearning(16:19);IntExp(8).deltanightFBbaseline=IntExp(6).deltanightFBlearning(14:15);
IntExp(8).FFmorningfitFBlearning=IntExp(6).FFmorningfitFBlearning(16:19);IntExp(8).FFmorningfitFBbaseline=IntExp(6).FFmorningfitFBlearning(14:15);
IntExp(8).AllDaysFBlearning=IntExp(6).AllDaysFBlearning(16:19);IntExp(8).AllDaysFBbaseline=IntExp(6).AllDaysFBlearning(14:15);
IntExp(8).FFeveningfitFBlearning=IntExp(6).FFeveningfitFBlearning(16:19);IntExp(8).FFeveningfitFBbaseline=IntExp(6).FFeveningfitFBlearning(14:15);
IntExp(8).deltadayCTlearning=IntExp(6).deltadayCTlearning(16:19);IntExp(8).deltadayCTbaseline=IntExp(6).deltadayCTlearning(14:15);
IntExp(8).deltanightCTlearning=IntExp(6).deltanightCTlearning(16:19);IntExp(8).deltanightCTbaseline=IntExp(6).deltanightCTlearning(14:15);
IntExp(8).FFmorningfitCTlearning=IntExp(6).FFmorningfitCTlearning(16:19);IntExp(8).FFmorningfitCTbaseline=IntExp(6).FFmorningfitCTlearning(14:15);
IntExp(8).AllDaysCTlearning=IntExp(6).AllDaysCTlearning(16:19);IntExp(8).AllDaysCTbaseline=IntExp(6).AllDaysCTlearning(14:15);
IntExp(8).FFeveningfitCTlearning=IntExp(6).FFeveningfitCTlearning(16:19);IntExp(8).FFeveningfitCTbaseline=IntExp(6).FFeveningfitCTlearning(14:15);


IntExp(6).deltadayFBlearning=IntExp(6).deltadayFBlearning(3:end);IntExp(6).deltadayFBbaseline=IntExp(6).deltadayFBlearning(1:2);
IntExp(6).deltanightFBlearning=IntExp(6).deltanightFBlearning(3:end);IntExp(6).deltanightFBbaseline=IntExp(6).deltanightFBlearning(1:2);
IntExp(6).FFmorningfitFBlearning=IntExp(6).FFmorningfitFBlearning(3:end);IntExp(6).FFmorningfitFBbaseline=IntExp(6).FFmorningfitFBlearning(1:2);
IntExp(6).AllDaysFBlearning=IntExp(6).AllDaysFBlearning(3:end);IntExp(6).AllDaysFBbaseline=IntExp(6).AllDaysFBlearning(1:2);
IntExp(6).FFeveningfitFBlearning=IntExp(6).FFeveningfitFBlearning(3:end);IntExp(6).FFeveningfitFBbaseline=IntExp(6).FFeveningfitFBlearning(1:2);
IntExp(6).deltadayCTlearning=IntExp(6).deltadayCTlearning(3:end);IntExp(6).deltadayCTbaseline=IntExp(6).deltadayCTlearning(1:2);
IntExp(6).deltanightCTlearning=IntExp(6).deltanightCTlearning(3:end);IntExp(6).deltanightCTbaseline=IntExp(6).deltanightCTlearning(1:2);
IntExp(6).FFmorningfitCTlearning=IntExp(6).FFmorningfitCTlearning(3:end);IntExp(6).FFmorningfitCTbaseline=IntExp(6).FFmorningfitCTlearning(1:2);
IntExp(6).AllDaysCTlearning=IntExp(6).AllDaysCTlearning(3:end);IntExp(6).AllDaysCTbaseline=IntExp(6).AllDaysCTlearning(1:2);
IntExp(6).FFeveningfitCTlearning=IntExp(6).FFeveningfitCTlearning(3:end);IntExp(6).FFeveningfitCTbaseline=IntExp(6).FFeveningfitCTlearning(1:2);

                    % TW data - already processed
                    %               %  load /bulbul1/LMAN_microstimulation/TWr87g80_sumdata.mat
                    %                 CTind=crct_com(1:53246);
                    %                 FBind=crfb_com(1:8919);
                    %                 vlsCT=vls(CTind,:);
                    %                 vlsFB=vls(FBind,:);
                    %                 vlsCTexp1=vlsCT(find(vlsCT(:,1)<734500),:);
                    %                 vlsFBexp1=vlsFB(find(vlsFB(:,1)<734500),:);
                    %                 vlsCTexp2=vlsCT(find(vlsCT(:,1)>734500),:);
                    %                 vlsFBexp2=vlsFB(find(vlsFB(:,1)>734500),:);
                    %                 [IntExp(5).deltadayFBlearning,IntExp(5).deltanightFBlearning,IntExp(5).AllDaysFBlearning,IntExp(5).FFmorningfitFBlearning,IntExp(5).FFeveningfitFBlearning]...
                    %                     =daynightTW(vlsFBexp1);
                    %                 [IntExp(6).deltadayFBlearning,IntExp(6).deltanightFBlearning,IntExp(6).AllDaysFBlearning,IntExp(6).FFmorningfitFBlearning,IntExp(6).FFeveningfitFBlearning]...
                    %                     =daynightTW(vlsFBexp2);
                    %                 [IntExp(5).deltadayCTlearning,IntExp(5).deltanightCTlearning,IntExp(5).AllDaysCTlearning,IntExp(5).FFmorningfitCTlearning,IntExp(5).FFeveningfitCTlearning]...
                    %                     =daynightTW(vlsCTexp1);
                    %                 [IntExp(6).deltadayCTlearning,IntExp(6).deltanightCTlearning,IntExp(6).AllDaysCTlearning,IntExp(6).FFmorningfitCTlearning,IntExp(6).FFeveningfitCTlearning]...
                    %                     =daynightTW(vlsCTexp2);
                    % 
                    %         IntExp(5).deltadayFBlearning=IntExp(5).deltadayFBlearning(4:end);IntExp(5).deltadayFBbaseline=IntExp(5).deltadayFBlearning(1:2); 
                    %         IntExp(5).deltanightFBlearning=IntExp(5).deltanightFBlearning(4:end);IntExp(5).deltanightFBbaseline=IntExp(5).deltanightFBlearning(1:2); 
                    %         IntExp(5).FFmorningfitFBlearning=IntExp(5).FFmorningfitFBlearning(4:end);IntExp(5).FFmorningfitFBbaseline=IntExp(5).FFmorningfitFBlearning(1:2); 
                    %         IntExp(5).AllDaysFBlearning=IntExp(5).AllDaysFBlearning(4:end);IntExp(5).AllDaysFBbaseline=IntExp(5).AllDaysFBlearning(1:2); 
                    %         IntExp(5).FFeveningfitFBlearning=IntExp(5).FFeveningfitFBlearning(4:end);IntExp(5).FFeveningfitFBbaseline=IntExp(5).FFeveningfitFBlearning(1:2); 
                    %         IntExp(5).deltadayCTlearning=IntExp(5).deltadayCTlearning(4:end);IntExp(5).deltadayCTbaseline=IntExp(5).deltadayCTlearning(1:2); 
                    %         IntExp(5).deltanightCTlearning=IntExp(5).deltanightCTlearning(4:end);IntExp(5).deltanightCTbaseline=IntExp(5).deltanightCTlearning(1:2); 
                    %         IntExp(5).FFmorningfitCTlearning=IntExp(5).FFmorningfitCTlearning(4:end);IntExp(5).FFmorningfitCTbaseline=IntExp(5).FFmorningfitCTlearning(1:2); 
                    %         IntExp(5).AllDaysCTlearning=IntExp(5).AllDaysCTlearning(4:end);IntExp(5).AllDaysCTbaseline=IntExp(5).AllDaysCTlearning(1:2); 
                    %         IntExp(5).FFeveningfitCTlearning=IntExp(5).FFeveningfitCTlearning(4:end);IntExp(5).FFeveningfitCTbaseline=IntExp(5).FFeveningfitCTlearning(1:2); 
                    %         IntExp(6).deltadayFBlearning=IntExp(6).deltadayFBlearning(3:end);IntExp(6).deltadayFBbaseline=IntExp(6).deltadayFBlearning(1:2); 
                    %         IntExp(6).deltanightFBlearning=IntExp(6).deltanightFBlearning(3:end);IntExp(6).deltanightFBbaseline=IntExp(6).deltanightFBlearning(1:2); 
                    %         IntExp(6).FFmorningfitFBlearning=IntExp(6).FFmorningfitFBlearning(3:end);IntExp(6).FFmorningfitFBbaseline=IntExp(6).FFmorningfitFBlearning(1:2); 
                    %         IntExp(6).AllDaysFBlearning=IntExp(6).AllDaysFBlearning(3:end);IntExp(6).AllDaysFBbaseline=IntExp(6).AllDaysFBlearning(1:2); 
                    %         IntExp(6).FFeveningfitFBlearning=IntExp(6).FFeveningfitFBlearning(3:end);IntExp(6).FFeveningfitFBbaseline=IntExp(6).FFeveningfitFBlearning(1:2); 
                    %         IntExp(6).deltadayCTlearning=IntExp(6).deltadayCTlearning(3:end);IntExp(6).deltadayCTbaseline=IntExp(6).deltadayCTlearning(1:2); 
                    %         IntExp(6).deltanightCTlearning=IntExp(6).deltanightCTlearning(3:end);IntExp(6).deltanightCTbaseline=IntExp(6).deltanightCTlearning(1:2); 
                    %         IntExp(6).FFmorningfitCTlearning=IntExp(6).FFmorningfitCTlearning(3:end);IntExp(6).FFmorningfitCTbaseline=IntExp(6).FFmorningfitCTlearning(1:2); 
                    %         IntExp(6).AllDaysCTlearning=IntExp(6).AllDaysCTlearning(3:end);IntExp(6).AllDaysCTbaseline=IntExp(6).AllDaysCTlearning(1:2); 
                    %         IntExp(6).FFeveningfitCTlearning=IntExp(6).FFeveningfitCTlearning(3:end);IntExp(6).FFeveningfitCTbaseline=IntExp(6).FFeveningfitCTlearning(1:2); 
                    % load /bulbul1/LMAN_microstimulation/NapStim/processedExps.mat
% 
%                         performance=ExpJC(6).performance;
%                         indFB=find(performance(:,2));
%                         indCT=find(~performance(:,2));
%                         clear vlsCTexp1 vlsFBexp1
%                         vlsCTexp1(:,1)=performance(indCT,1);
%                         vlsCTexp1(:,2)=2*performance(indCT,3);
%                         vlsFBexp1(:,1)=performance(indFB,1);
%                         vlsFBexp1(:,2)=2*performance(indFB,3);
%                         [IntExp(7).deltadayFBlearning,IntExp(7).deltanightFBlearning,IntExp(7).AllDaysFBlearning,IntExp(7).FFmorningfitFBlearning,IntExp(7).FFeveningfitFBlearning]...
%                             =daynightTW(vlsFBexp1);
%                         [IntExp(7).deltadayCTlearning,IntExp(7).deltanightCTlearning,IntExp(7).AllDaysCTlearning,IntExp(7).FFmorningfitCTlearning,IntExp(7).FFeveningfitCTlearning]...
%                             =daynightTW(vlsCTexp1);
%                         %%%
%                         performance=ExpJC(8).performance;
%                         indFB=find(performance(:,2));
%                         indCT=find(~performance(:,2));
%                         clear vlsCTexp1 vlsFBexp1
%                         vlsCTexp1(:,1)=performance(indCT,1);
%                         vlsCTexp1(:,2)=2*performance(indCT,3);
%                         vlsFBexp1(:,1)=performance(indFB,1);
%                         vlsFBexp1(:,2)=2*performance(indFB,3);
%                         [IntExp(8).deltadayFBlearning,IntExp(8).deltanightFBlearning,IntExp(8).AllDaysFBlearning,IntExp(8).FFmorningfitFBlearning,IntExp(8).FFeveningfitFBlearning]...
%                             =daynightTW(vlsFBexp1);
%                         [IntExp(8).deltadayCTlearning,IntExp(8).deltanightCTlearning,IntExp(8).AllDaysCTlearning,IntExp(8).FFmorningfitCTlearning,IntExp(8).FFeveningfitCTlearning]...
%                             =daynightTW(vlsCTexp1);
% 

clear all
load /bulbul1/IntExp030612.mat

%%%%%%%%%%%%%%%%%%%%%
%%% Average trajectories
%%%%%%%%%%%%%%%%%%%%%
            clear TrajST TrajCT mnCT mnST seCT seST seA
            for i=[2:5 7:8]
                for k=1:length(IntExp(i).deltanightCTlearning)
                    TrajCT(i-1,k*2-1)=IntExp(i).dir*(IntExp(i).deltadayCTlearning(k));
                    TrajCT(i-1,k*2)=IntExp(i).dir*IntExp(i).deltanightCTlearning(k);
                end
                for k=1:length(IntExp(i).deltanightCTlearning)
                    TrajST(i-1,k*2-1)=IntExp(i).dir*(IntExp(i).deltadayFBlearning(k));
                    TrajST(i-1,k*2)=IntExp(i).dir*IntExp(i).deltanightFBlearning(k);
                end
            end
            cTrajC=cumsum(TrajCT');
            cTrajS=cumsum(TrajST');
            mnCT(1:4)=mean(cTrajC(1:4,[1:4 6:7])');mnST(1:4)=mean(cTrajS(1:4,[1:4 6:7])');mnCT(5:6)=mean(cTrajC(5:6,[1:4 7])');mnST(5:6)=mean(cTrajS(5:6,[1:4 7])');mnCT(7:8)=mean(cTrajC(7:8,[1:4 ])');mnST(7:8)=mean(cTrajS(7:8,[1:4 ])');
            seCT(1:4)=std(cTrajC(1:4,[1:4 6:7])')/sqrt(6);seST(1:4)=std(cTrajS(1:4,[1:4 6:7])')/sqrt(6);seCT(5:6)=std(cTrajC(5:6,[1:4 7])')/sqrt(5);seST(5:6)=std(cTrajS(5:6,[1:4 7])')/sqrt(5);seCT(7:8)=std(cTrajC(7:8,[1:4 ])')/sqrt(4);seST(7:8)=std(cTrajS(7:8,[1:4 ])')/sqrt(4);
            seA(1:4)=std([cTrajC(1:4,[1:4 6:7])-cTrajS(1:4,[1:4 6:7])]')/sqrt(6);seA(5:6)=std([cTrajC(5:6,[1:4 7])-cTrajS(5:6,[1:4 7])]')/sqrt(5);seA(7:8)=std([cTrajC(7:8,[1:4 ])-cTrajS(7:8,[1:4 ])]')/sqrt(4);
            
            figure;hold on;
            plot([0 mnCT],'Linewidth',2)
            plot([0 mnST],'r','Linewidth',2)
            plot([2 3],mnCT(1:2),'k','Linewidth',2);plot([2 3],mnST(1:2),'k','Linewidth',2)
            plot([4 5],mnCT(3:4),'k','Linewidth',2);plot([4 5],mnST(3:4),'k','Linewidth',2)
            plot([6 7],mnCT(5:6),'k','Linewidth',2);plot([6 7],mnST(5:6),'k','Linewidth',2)
            
            mnA=mnCT-mnST;
            plot([0 mnA],'g')
            plot([2 3],mnA(1:2),'k','Linewidth',2);plot([4 5],mnA(3:4),'k','Linewidth',2);plot([6 7],mnA(5:6),'k','Linewidth',2)
            
            for k=1:8
                plot([k+1 k+1],[mnCT(k)+seCT(k);mnCT(k)-seCT(k)],'b')
                plot([k+1 k+1],[mnST(k)+seST(k);mnST(k)-seST(k)],'r')
                plot([k+1 k+1],[mnA(k)+seA(k);mnA(k)-seA(k)],'g')
            end
            xlim([1 8])
%%%% Individual experiments
figure;hold on;subplot(122);hold on;
    clr='rgbyamc';
    for j=[1:4 7]
    plot([0 cTrajS(1:8,j)'],'Linewidth',2,'color',clr(j))
    plot([2 3],cTrajS(1:2,j),'k','Linewidth',2);plot([2 3],cTrajS(1:2,j),'k','Linewidth',2)
    plot([4 5],cTrajS(3:4,j),'k','Linewidth',2);plot([4 5],cTrajS(3:4,j),'k','Linewidth',2)
    plot([6 7],cTrajS(5:6,j),'k','Linewidth',2);plot([6 7],cTrajS(5:6,j),'k','Linewidth',2)
    end
    j=6;
    plot([0 cTrajS(1:4,j)'],'Linewidth',2,'color',clr(j))
    plot([2 3],cTrajS(1:2,j),'k','Linewidth',2);plot([2 3],cTrajS(1:2,j),'k','Linewidth',2)
    plot([4 5],cTrajS(3:4,j),'k','Linewidth',2);plot([4 5],cTrajS(3:4,j),'k','Linewidth',2)
    xlim([1 8]);ylim([-10 140])
subplot(121);hold on;
    for j=[1:4 7]
    plot([0 cTrajC(1:8,j)'],'Linewidth',2,'color',clr(j))
    plot([2 3],cTrajC(1:2,j),'k','Linewidth',2);plot([2 3],cTrajC(1:2,j),'k','Linewidth',2)
    plot([4 5],cTrajC(3:4,j),'k','Linewidth',2);plot([4 5],cTrajC(3:4,j),'k','Linewidth',2)
    plot([6 7],cTrajC(5:6,j),'k','Linewidth',2);plot([6 7],cTrajC(5:6,j),'k','Linewidth',2)
    end
    j=6;
    plot([0 cTrajC(1:4,j)'],'Linewidth',2,'color',clr(j))
    plot([2 3],cTrajC(1:2,j),'k','Linewidth',2);plot([2 3],cTrajC(1:2,j),'k','Linewidth',2)
    plot([4 5],cTrajC(3:4,j),'k','Linewidth',2);plot([4 5],cTrajC(3:4,j),'k','Linewidth',2)
    xlim([1 8]);ylim([-10 140])



%%%%%%%%%%%%%%%%%%%%%
%%% Bar plots - all
%%%%%%%%%%%%%%%%%%%%%


  % Only consider the first four days
  	lastday=[4 4 3 4 4 3 2 3];
  % For AFP-dep, only consider the first two days
    indfirst2=[1:2 5:6 8:9 12:13 16:17 18:19];
    indfirst2=[1 5 8 12 16 18];
    
  %%%%%%%%
            %  STIM
                        clear LearnST 
                        for i=1:length(lastday)
                            % 1=all, 2=day, 3=night
                                LearnST{i}{1}=IntExp(i).dir*(IntExp(i).deltadayFBlearning(1:lastday(i))+IntExp(i).deltanightFBlearning(1:lastday(i)));
                                LearnST{i}{2}=IntExp(i).dir*IntExp(i).deltadayFBlearning(1:lastday(i));
                                LearnST{i}{3}=IntExp(i).dir*IntExp(i).deltanightFBlearning(1:lastday(i));    
                        end

                        allST=([LearnST{2}{1} LearnST{3}{1} LearnST{4}{1} LearnST{5}{1} LearnST{7}{1} LearnST{8}{1}]);
                        dayST=([LearnST{2}{2} LearnST{3}{2} LearnST{4}{2} LearnST{5}{2} LearnST{7}{2} LearnST{8}{2}]);% LearnST{6}{2}]);
                        nightST=([LearnST{2}{3} LearnST{3}{3} LearnST{4}{3} LearnST{5}{3} LearnST{7}{3} LearnST{8}{3}]);% LearnST{6}{3}]);
            %   CATCH
                        clear LearnCT 
                        for i=1:length(lastday)
                                LearnCT{i}{1}=IntExp(i).dir*(IntExp(i).deltadayCTlearning(1:lastday(i))+IntExp(i).deltanightCTlearning(1:lastday(i)));
                                LearnCT{i}{2}=IntExp(i).dir*IntExp(i).deltadayCTlearning(1:lastday(i));
                                LearnCT{i}{3}=IntExp(i).dir*IntExp(i).deltanightCTlearning(1:lastday(i));    
                        end
                        %
                        allCT=([LearnCT{2}{1} LearnCT{3}{1} LearnCT{4}{1} LearnCT{5}{1} LearnCT{7}{1} LearnCT{8}{1} ]);% LearnCT{6}{1}]);
                        dayCT=([LearnCT{2}{2} LearnCT{3}{2} LearnCT{4}{2} LearnCT{5}{2} LearnCT{7}{2} LearnCT{8}{2}]);% LearnCT{6}{2}]);
                        nightCT=([LearnCT{2}{3} LearnCT{3}{3} LearnCT{4}{3} LearnCT{5}{3} LearnCT{7}{3} LearnCT{8}{3}]);% LearnCT{6}{3}]);
            % DAYS (only considering the first four days, i.e. before asymptote)
                    figure;hold on;
                    % total
                        bar(1,mean(dayCT));bar(2,mean(nightCT))
                        plot([1 1],[mean(dayCT)-std(dayCT)/sqrt(length(dayCT)) mean(dayCT)+std(dayCT)/sqrt(length(dayCT))],'-');
                        plot([2 2],[mean(nightCT)-std(nightCT)/sqrt(length(nightCT)) mean(nightCT)+std(nightCT)/sqrt(length(nightCT))],'-');
                        plot(1,dayCT,'o')
                        plot(2,nightCT,'o')
                    % motor pathway
                        bar(4,mean(dayST));bar(5,mean(nightST))
                        plot([4 4],[mean(dayST)-std(dayST)/sqrt(length(dayST)) mean(dayST)+std(dayST)/sqrt(length(dayST))],'-');
                        plot([5 5],[mean(nightST)-std(nightST)/sqrt(length(nightST)) mean(nightST)+std(nightST)/sqrt(length(nightST))],'-');
                        plot(4,dayST,'o')
                        plot(5,nightST,'o')

                    % AFP
                        bar(7,mean(dayCT)-mean(dayST));bar(8,mean(nightCT)-mean(nightST))
                        plot([7 7],[mean(dayCT-dayST)-std(dayCT-dayST)/sqrt(length(dayCT-dayST)) mean(dayCT-dayST)+std(dayCT-dayST)/sqrt(length(dayCT-dayST))],'-');
                        plot([8 8],[mean(nightCT-nightST)-std(nightCT-nightST)/sqrt(length(nightCT-nightST)) mean(nightCT-nightST)+std(nightCT-nightST)/sqrt(length(nightCT-nightST))],'-');
                    % AFP- just the first two days
                        bar(10,mean(dayCT(indfirst2))-mean(dayST(indfirst2)));bar(11,mean(nightCT(indfirst2))-mean(nightST(indfirst2)))
                        plot([10 10],[mean(dayCT(indfirst2)-dayST(indfirst2))-std(dayCT(indfirst2)-dayST(indfirst2))/sqrt(length(dayCT(indfirst2)-dayST(indfirst2))) mean(dayCT(indfirst2)-dayST(indfirst2))+std(dayCT(indfirst2)-dayST(indfirst2))/sqrt(length(dayCT(indfirst2)-dayST(indfirst2)))],'-');
                        plot([11 11],[mean(nightCT(indfirst2)-nightST(indfirst2))-std(nightCT(indfirst2)-nightST(indfirst2))/sqrt(length(nightCT(indfirst2)-nightST(indfirst2))) mean(nightCT(indfirst2)-nightST(indfirst2))+std(nightCT(indfirst2)-nightST(indfirst2))/sqrt(length(nightCT(indfirst2)-nightST(indfirst2)))],'-');
            
                    
%%%%%%%%%%%%%%%%%%                    
% What do things look like first thing in the morning?
%%%%%%%%%%%%%%%%%%
                clear hoursCT dir daydeltCT nightdeltCT
                count=0;                
                for i=[2:5 7:8]
                    for k=1:lastday(i);
                        count=count+1;
                        dir(count)=IntExp(i).dir;
                        day(count)=k;
                        %
                        if length(IntExp(i).AllDaysCTlearning)>k
                            daydeltCT(count)=dir(count)*(mean(IntExp(i).AllDaysCTlearning(k).FFdata(end-50:end))-mean(IntExp(i).AllDaysCTlearning(k).FFdata(1:50)));
                            nightdeltCT(count)=dir(count)*(mean(IntExp(i).AllDaysCTlearning(k+1).FFdata(1:50))-mean(IntExp(i).AllDaysCTlearning(k).FFdata(end-50:end)));
                        else
                            daydeltCT(count)=nan;
                            nightdeltCT(count)=nan;
                        end
                        %
                        if i<5
                            tv=(IntExp(i).AllDaysCTlearning(k).tdata-min(IntExp(i).AllDaysCTlearning(k).tdata))/3600;
                        else
                            tv=(IntExp(i).AllDaysCTlearning(k).tdata-min(IntExp(i).AllDaysCTlearning(k).tdata))*24;
                        end
                        for h=1:14
                            hoursCT(h,count)=mean(IntExp(i).AllDaysCTlearning(k).FFdata(find(tv>h-1 & tv<h)));
                        end
                        if length(IntExp(i).AllDaysCTlearning)>k
                            if i<5
                                tv=(IntExp(i).AllDaysCTlearning(k+1).tdata-min(IntExp(i).AllDaysCTlearning(k+1).tdata))/3600;
                            else
                                tv=(IntExp(i).AllDaysCTlearning(k+1).tdata-min(IntExp(i).AllDaysCTlearning(k+1).tdata))*24;
                            end
                            hoursCT(15,count)=mean(IntExp(i).AllDaysCTlearning(k+1).FFdata(find(tv>0 & tv<1)));
                        else
                            hoursCT(15,count)=NaN;
                        end
                    end    
                end
                
                for i=1:size(hoursCT,2)
                    hoursCT(:,i)=dir(i)*(hoursCT(:,i)-hoursCT(1,i));
                end
                
                clear mnhrsCT sehrsCT
                for i=1:size(hoursCT,1)
                    mnhrsCT(i)=mean(hoursCT(i,~isnan(hoursCT(i,:))));
                    sehrsCT(i)=std(hoursCT(i,~isnan(hoursCT(i,:))))/sqrt(length(~isnan(hoursCT(i,:))));
                end
               
                clear hoursST daydeltST nightdeltST
                count=0;                
                for i=[2:5 7:8]
                    for k=1:lastday(i);
                        count=count+1;
                       dir(count)=IntExp(i).dir;
                        %
                        %
                        if length(IntExp(i).AllDaysFBlearning)>k
                            daydeltST(count)=dir(count)*(mean(IntExp(i).AllDaysFBlearning(k).FFdata(end-10:end))-mean(IntExp(i).AllDaysFBlearning(k).FFdata(1:10)));
                            nightdeltST(count)=dir(count)*(mean(IntExp(i).AllDaysFBlearning(k+1).FFdata(1:10))-mean(IntExp(i).AllDaysFBlearning(k).FFdata(end-10:end)));
                        else
                            daydeltST(count)=nan;
                            nightdeltST(count)=nan;
                        end
                        %
                        %
                        if i<5
                            tv=(IntExp(i).AllDaysFBlearning(k).tdata-min(IntExp(i).AllDaysFBlearning(k).tdata))/3600;
                        else
                            tv=(IntExp(i).AllDaysFBlearning(k).tdata-min(IntExp(i).AllDaysFBlearning(k).tdata))*24;
                        end
                        for h=1:14
                            hoursST(h,count)=mean(IntExp(i).AllDaysFBlearning(k).FFdata(find(tv>h-1 & tv<h)));
                        end
                        if length(IntExp(i).AllDaysCTlearning)>k
                            if i<5
                                tv=(IntExp(i).AllDaysFBlearning(k+1).tdata-min(IntExp(i).AllDaysFBlearning(k+1).tdata))/3600;
                            else
                                tv=(IntExp(i).AllDaysFBlearning(k+1).tdata-min(IntExp(i).AllDaysFBlearning(k+1).tdata))*24;
                            end
                            hoursST(15,count)=mean(IntExp(i).AllDaysFBlearning(k+1).FFdata(find(tv>0 & tv<1)));
                        else
                            hoursST(15,count)=NaN;
                        end

                    end 
                end

                for i=1:size(hoursST,2)
                    hoursST(:,i)=dir(i)*(hoursST(:,i)-hoursST(1,i));
                end
                
                clear mnhrsST sehrsST
                for i=1:size(hoursST,1)
                    mnhrsST(i)=mean(hoursST(i,~isnan(hoursST(i,:))));
                    sehrsST(i)=std(hoursST(i,~isnan(hoursST(i,:))))/sqrt(length(~isnan(hoursST(i,:))));
                end
           % Define sleep changes 
                mean(hoursST(14,[1 3 6 8:15 18 19 20])-hoursST(15,[1 3 6 8:15 18 19 20]))
                std(hoursST(14,[1 3 6 8:15 18 19 20])-hoursST(15,[1 3 6 8:15 18 19 20]))/sqrt(14)
                
           %%% hourly     
                figure;hold on;
                plot(mnhrsCT([1:14]),'c');plot(mnhrsST([1:14]),'r')
                plot(18,mnhrsCT(15),'co');plot(18,mnhrsST(15),'ro')
                
                for k=1:14
                    plot([k k],[mnhrsCT(k)+sehrsCT(k);mnhrsCT(k)-sehrsCT(k)],'c')
                    plot([k k],[mnhrsST(k)+sehrsST(k);mnhrsST(k)-sehrsST(k)],'r')
                end
                k=15;
                    plot([k+3 k+3],[mnhrsCT(k)+sehrsCT(k);mnhrsCT(k)-sehrsCT(k)],'c')
                    plot([k+3 k+3],[mnhrsST(k)+sehrsST(k);mnhrsST(k)-sehrsST(k)],'r')
                plot([14 18],mnhrsCT([14 15]),'c');plot([14 18],mnhrsST([14 15]),'r')
                plot([0 19],[0 0],'k')
                
                
                
                
 %%%%
 % correlations
 %%%%
 
    % end of day (1hrs) to morning (1hrs)
            %         for i=1:size(hoursCT,2)   % for each experiment
            %             reals=hoursCT(find(~isnan(hoursCT(1:14,i))),i);
            %             deltCTday(i)=mean([reals(end) ]);
            %             deltCTnight(i)=hoursCT(15,i)-deltCTday(i);
            %         end
            %         for i=1:size(hoursST,2)   % for each experiment
            %             reals=hoursST(find(~isnan(hoursST(1:14,i))),i);
            %             deltSTday(i)=mean([reals(end)]);
            %             deltSTnight(i)=hoursST(15,i)-deltCTday(i);
            %         end
            %         figure;hold on;subplot(1,2,1);hold on;plot(deltCTday,deltCTnight,'o');plot([-60 60],[0 0]);plot([0 0],[-60 60]);xlim([-60 60]);ylim([-60 60])                
            %         subplot(122);hold on;plot(deltSTday,deltSTnight,'ro');plot([-60 60],[0 0]);plot([0 0],[-60 60]);xlim([-60 60]);ylim([-60 60])                
   %%% % try it by number of songs instead
         figure;hold on;subplot(1,2,1);hold on;plot(daydeltCT,nightdeltCT,'o');plot([-100 100],[0 0]);plot([0 0],[-100 100]);xlim([-100 100]);ylim([-100 100])                
        subplot(122);hold on;plot(daydeltST,nightdeltST,'ro');plot([-100 100],[0 0]);plot([0 0],[-100 100]);xlim([-100 100]);ylim([-100 100])                
   
        
        
        
        
        
%%%% bar graphs

 %%% end of day (1hrs) to morning (1hrs)
                % clear deltCTday deltCTnight deltSTday deltSTnight
                %          for i=1:size(hoursCT,2)   % for each experiment
                %             reals=hoursCT(find(~isnan(hoursCT(1:14,i))),i);
                %             deltCTday(i)=mean([reals(end)]);
                %             deltCTnight(i)=hoursCT(15,i)-deltCTday(i);
                %         end
                %         for i=1:size(hoursST,2)   % for each experiment
                %             reals=hoursST(find(~isnan(hoursST(1:14,i))),i);
                %             deltSTday(i)=mean([reals(end) ]);
                %             deltSTnight(i)=hoursST(15,i)-deltCTday(i);
                %         end
                %         figure;hold on;
                %         bar(1,mean(deltCTday));bar(2,mean(deltCTnight(~isnan(deltSTnight))));
                %         bar(4,mean(deltSTday));bar(5,mean(deltSTnight(~isnan(deltSTnight))))
                %         plot([1 1],[mean(deltCTday)-std(deltCTday)/sqrt(length(deltCTday)) mean(deltCTday)+std(deltCTday)/sqrt(length(deltCTday))])
                %         plot([2 2],[mean(deltCTnight(~isnan(deltSTnight)))-std(deltCTnight(~isnan(deltSTnight)))/sqrt(length(deltCTnight(~isnan(deltSTnight)))) mean(deltCTnight(~isnan(deltSTnight)))+std(deltCTnight(~isnan(deltSTnight)))/sqrt(length(deltCTnight(~isnan(deltSTnight))))])
                %         plot([4 4],[mean(deltSTday)-std(deltSTday)/sqrt(length(deltSTday)) mean(deltSTday)+std(deltSTday)/sqrt(length(deltSTday))])
                %         plot([5 5],[mean(deltSTnight(~isnan(deltSTnight)))-std(deltSTnight(~isnan(deltSTnight)))/sqrt(length(deltSTnight(~isnan(deltSTnight)))) mean(deltSTnight(~isnan(deltSTnight)))+std(deltSTnight(~isnan(deltSTnight)))/sqrt(length(deltSTnight(~isnan(deltSTnight))))])
%%% try it by number of songs instead
        figure;hold on;
        bar(1,mean(daydeltCT(~isnan(daydeltCT))));bar(2,mean(nightdeltCT(~isnan(nightdeltST))));
        bar(4,mean(daydeltST(~isnan(daydeltST))));bar(5,mean(nightdeltST(~isnan(nightdeltST))))
        plot([1 1],[mean(daydeltCT(~isnan(daydeltCT)))-std(daydeltCT(~isnan(daydeltCT)))/sqrt(length(daydeltCT(~isnan(daydeltCT)))) mean(daydeltCT(~isnan(daydeltCT)))+std(daydeltCT(~isnan(daydeltCT)))/sqrt(length(daydeltCT(~isnan(daydeltCT))))])
        plot([2 2],[mean(nightdeltCT(~isnan(nightdeltST)))-std(nightdeltCT(~isnan(nightdeltST)))/sqrt(length(nightdeltCT(~isnan(nightdeltST)))) mean(nightdeltCT(~isnan(nightdeltST)))+std(nightdeltCT(~isnan(nightdeltST)))/sqrt(length(nightdeltCT(~isnan(nightdeltST))))])
        plot([4 4],[mean(daydeltST(~isnan(daydeltST)))-std(daydeltST(~isnan(daydeltST)))/sqrt(length(daydeltST(~isnan(daydeltST)))) mean(daydeltST(~isnan(daydeltST)))+std(daydeltST(~isnan(daydeltST)))/sqrt(length(daydeltST(~isnan(daydeltST))))])
        plot([5 5],[mean(nightdeltST(~isnan(nightdeltST)))-std(nightdeltST(~isnan(nightdeltST)))/sqrt(length(nightdeltST(~isnan(nightdeltST)))) mean(nightdeltST(~isnan(nightdeltST)))+std(nightdeltST(~isnan(nightdeltST)))/sqrt(length(nightdeltST(~isnan(nightdeltST))))])





%%%%
 % daily
%%%%   

%%%%%%%%%%%%%%%%%%%%%
%%% Average trajectories
%%%%%%%%%%%%%%%%%%%%%


%%% hours
            %             clear TrajSTh TrajCTh mnCT mnST 
            %             for i=1:4
            %                 TrajCTh(i*2-1,1:length(find(day==i)))=deltCTday(day==i);
            %                 TrajSTh(i*2-1,1:length(find(day==i)))=deltSTday(day==i);
            %                 TrajCTh(i*2,1:length(find(day==i)))=deltCTnight(day==i);
            %                 TrajSTh(i*2,1:length(find(day==i)))=deltSTnight(day==i);
            %             end
            %             
            %             cTrajCh=cumsum(TrajCTh);
            %             cTrajSh=cumsum(TrajSTh);
            %             mnCT(1:3)=mean(cTrajCh(1:3,:)');mnCT(4)=mean(cTrajCh(4,[1:4 6])');mnCT(5:6)=mean(cTrajCh(5:6,[1:4])');mnCT(7:8)=mean(cTrajCh(7:8,[1:3])');
            %             mnST(1:3)=mean(cTrajSh(1:3,:)');mnST(4)=mean(cTrajSh(4,[1:4 6])');mnST(5:6)=mean(cTrajSh(5:6,[1:4])');mnST(7:8)=mean(cTrajSh(7:8,[1:3])');
            %             figure;hold on;plot(mnCT);plot(mnST,'r')
            %          
            %            seCT(1:3)=std(cTrajCh(1:3,:)')/sqrt(6);seCT(4)=std(cTrajCh(4,[1:4 6])')/sqrt(5);seCT(5:6)=std(cTrajCh(5:6,[1:4])')/sqrt(4);seCT(7:8)=std(cTrajCh(7:8,[1:3])')/sqrt(3);
            %             seST(1:3)=std(cTrajSh(1:3,:)')/sqrt(6);seST(4)=std(cTrajSh(4,[1:4 6])')/sqrt(5);seST(5:6)=std(cTrajSh(5:6,[1:4])')/sqrt(4);seST(7:8)=std(cTrajSh(7:8,[1:3])')/sqrt(3);
            % 
            %             figure;hold on;
            %             plot([0 mnCT],'Linewidth',2)
            %             plot([0 mnST],'r','Linewidth',2)
            %             plot([2 3],mnCT(1:2),'k','Linewidth',2);plot([2 3],mnST(1:2),'k','Linewidth',2)
            %             plot([4 5],mnCT(3:4),'k','Linewidth',2);plot([4 5],mnST(3:4),'k','Linewidth',2)
            %             plot([6 7],mnCT(5:6),'k','Linewidth',2);plot([6 7],mnST(5:6),'k','Linewidth',2)
            %                         
            %             for k=1:8
            %                 plot([k+1 k+1],[mnCT(k)+seCT(k);mnCT(k)-seCT(k)],'b')
            %                 plot([k+1 k+1],[mnST(k)+seST(k);mnST(k)-seST(k)],'r')
            %             end
            %             xlim([1 8])
% try it by number of songs instead
             clear TrajSTh TrajCTh mnCT mnST cTrajCh cTrajSh
            for i=1:4
                TrajCTh(i*2-1,1:length(find(day==i)))=daydeltCT(day==i);
                TrajSTh(i*2-1,1:length(find(day==i)))=daydeltST(day==i);
                TrajCTh(i*2,1:length(find(day==i)))=nightdeltCT(day==i);
                TrajSTh(i*2,1:length(find(day==i)))=nightdeltST(day==i);
            end
            
            cTrajCh=cumsum(TrajCTh);
            cTrajSh=cumsum(TrajSTh);
            mnCT(1:2)=mean(cTrajCh(1:2,:)');mnCT(3:4)=mean(cTrajCh(4,[1:4 6])');mnCT(5:6)=mean(cTrajCh(5:6,[1:4])');mnCT(7:8)=mean(cTrajCh(7:8,[1:3])');
            mnST(1:2)=mean(cTrajSh(1:2,:)');mnST(3:4)=mean(cTrajSh(4,[1:4 6])');mnST(5:6)=mean(cTrajSh(5:6,[1:4])');mnST(7:8)=mean(cTrajSh(7:8,[1:3])');
            figure;hold on;plot(mnCT);plot(mnST,'r')
         
           seCT(1:2)=std(cTrajCh(1:2,:)')/sqrt(6);seCT(3:4)=std(cTrajCh(4,[1:4 6])')/sqrt(5);seCT(5:6)=std(cTrajCh(5:6,[1:4])')/sqrt(4);seCT(7:8)=std(cTrajCh(7:8,[1:3])')/sqrt(3);
            seST(1:2)=std(cTrajSh(1:2,:)')/sqrt(6);seST(3:4)=std(cTrajSh(4,[1:4 6])')/sqrt(5);seST(5:6)=std(cTrajSh(5:6,[1:4])')/sqrt(4);seST(7:8)=std(cTrajSh(7:8,[1:3])')/sqrt(3);

            figure;hold on;
            plot([0 mnCT],'Linewidth',2)
            plot([0 mnST],'r','Linewidth',2)
            plot([2 3],mnCT(1:2),'k','Linewidth',2);plot([2 3],mnST(1:2),'k','Linewidth',2)
            plot([4 5],mnCT(3:4),'k','Linewidth',2);plot([4 5],mnST(3:4),'k','Linewidth',2)
            plot([6 7],mnCT(5:6),'k','Linewidth',2);plot([6 7],mnST(5:6),'k','Linewidth',2)
                        
            for k=1:8
                plot([k+1 k+1],[mnCT(k)+seCT(k);mnCT(k)-seCT(k)],'b')
                plot([k+1 k+1],[mnST(k)+seST(k);mnST(k)-seST(k)],'r')
            end
            xlim([1 8])

       
        
        
        
        
 %%%%
 % example
 %%%%
             load /bulbul1/LMAN_microstimulation/Data020912.mat

              ExpIndex=4;
            figure;hold on;
            catchfrac1=5;catchfrac2=1;ravgwin=50;
            bct=[];bst=[];fct=[];fst=[];
            for i=36:45%length(experiment)
                window=[400:440];
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                bct=[bct b];
                timedata=b;
                FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                fct=[fct FFdata];
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if endpoints(j)>startpoints(j)+catchfrac-1
                        inds=[startpoints(j):catchfrac1:endpoints(j)];
                        plot(timedata(inds),FFdata(inds),'.','markersize',7,'Color','b')
                    end
                end
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                bst=[bst b];
                timedata=b;
                FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                fst=[fst FFdata];
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if endpoints(j)>startpoints(j)+catchfrac-1
                        inds=[startpoints(j):catchfrac2:endpoints(j)];           
                        plot(timedata(inds),FFdata(inds),'.','markersize',7,'Color','r')
                    end
                end
            end
            
     %%% HOURS
            lastpoints=find(diff(bct)>10);
            startpoints=[1 lastpoints+1];
            clear night dawn
        % by hour
%             for i=1:length(lastpoints)
%                 night(i)=mean(fct(bct>bct(lastpoints(i))-1 & bct<bct(lastpoints(i))));
%                 dawn(i)=mean(fct(bct>bct(startpoints(i)) & bct<bct(startpoints(i))+1));
%             end
        % by number of songs performed
            for i=1:length(lastpoints)
                night(i)=mean(fct(lastpoints(i)-50:lastpoints(i)));
                dawn(i)=mean(fct(startpoints(i):startpoints(i)+50));
            end

            for i=1:length(night)
                plot([bct(startpoints(i));bct(lastpoints(i))],[dawn(i);night(i)],'b-','Linewidth',4)
            end
         
            lastpoints=find(diff(bst)>10);
            startpoints=[1 lastpoints+1];
            clear night dawn
         % by hour
%             for i=1:length(lastpoints)
%                 night(i)=mean(fst(bst>bst(lastpoints(i))-1 & bst<bst(lastpoints(i))));
%                 dawn(i)=mean(fst(bst>bst(startpoints(i)) & bst<bst(startpoints(i))+1));
%             end
         % by number of songs performed
            for i=1:length(lastpoints)
                night(i)=mean(fct(lastpoints(i)-10:lastpoints(i)));
                dawn(i)=mean(fct(startpoints(i):startpoints(i)+10));
            end

            for i=1:length(night)
                plot([bst(startpoints(i));bst(lastpoints(i))],[dawn(i);night(i)],'r-','Linewidth',4)
            end
            
            dayST=106.1;
            nightST=38.3;
          
                
                
                
%%%%%%%%%%%%%%%%%%                    
% When was the hit rate adjusted
%%%%%%%%%%%%%%%%%%
                clear hoursCT dir hoursHitRate
                count=0;                
                for i=[2:5 7:8]
                    for k=1:lastday(i);
                        count=count+1;
                        dir(count)=IntExp(i).dir;
                        if i<5
                            tv=(IntExp(i).AllDaysCTlearning(k).tdata-min(IntExp(i).AllDaysCTlearning(k).tdata))/3600;
                            for h=1:14
                                fvmatrix=IntExp(i).AllDaysCTlearning(k).fv(find(tv>h-1 & tv<h));
                                if length(fvmatrix)>0
                                count2=0;
                                clear istrig
                                for ii=1:length(fvmatrix);
                                        count2=count2+1;
                                        istrig(ii)=fvmatrix(ii).TRIG(1);
                                end
                                hoursHitRate(h,count)=sum(istrig)/length(istrig);
                                end
                            end
                        else
                            tv=(IntExp(i).AllDaysCTlearning(k).tdata-min(IntExp(i).AllDaysCTlearning(k).tdata))*24;
                            for h=1:14
                                hrind=find(tv>h-1 & tv<h);
                                if length(hrind)>0
                                    hoursHitRate(h,count)=sum(IntExp(i).AllDaysCTlearning(k).hrdata(hrind))/length(hrind);
                                end
                            end
                        end
                    end    
                end
                
                clear mnhrsHR
                for i=1:size(hoursHitRate,1)
                    mnhrsHR(i)=mean(hoursHitRate(i,~isnan(hoursHitRate(i,:))));
                    sehrsHR(i)=std(hoursHitRate(i,~isnan(hoursHitRate(i,:))))/sqrt(length(~isnan(hoursHitRate(i,:))));
                end

           
            figure;hold on;
            plot(mnhrsHR)
            plot(mnhrsHR-sehrsHR)
            plot(mnhrsHR+sehrsHR)
            ylim([0 1])
            
            
            
            
            
            
    % BIRDS
        
%  clear allCT allST dayCT dayST nightCT nightST
% for k=2:6
%     allCT(k-1)=mean(LearnCT{k}{1});
%     allST(k-1)=mean(LearnST{k}{1});
%     dayCT(k-1)=mean(LearnCT{k}{2});
%     dayST(k-1)=mean(LearnST{k}{2});
%     nightCT(k-1)=mean(LearnCT{k}{3});
%     nightST(k-1)=mean(LearnST{k}{3});
% end
% 
%             
%         figure;hold on;
%         % total
%             bar(1,mean(dayCT));bar(2,mean(nightCT))
%             plot([1 1],[mean(dayCT)-std(dayCT)/sqrt(length(dayCT)) mean(dayCT)+std(dayCT)/sqrt(length(dayCT))],'-');
%             plot([2 2],[mean(nightCT)-std(nightCT)/sqrt(length(nightCT)) mean(nightCT)+std(nightCT)/sqrt(length(nightCT))],'-');
% 
%         % motor pathway
%             bar(4,mean(dayST));bar(5,mean(nightST))
%             plot([4 4],[mean(dayST)-std(dayST)/sqrt(length(dayST)) mean(dayST)+std(dayST)/sqrt(length(dayST))],'-');
%             plot([5 5],[mean(nightST)-std(nightST)/sqrt(length(nightST)) mean(nightST)+std(nightST)/sqrt(length(nightST))],'-');
% 
%         % AFP
%             bar(7,mean(dayCT)-mean(dayST));bar(8,mean(nightCT)-mean(nightST))
%             plot([7 7],[mean(dayCT-dayST)-std(dayCT-dayST)/sqrt(length(dayCT-dayST)) mean(dayCT-dayST)+std(dayCT-dayST)/sqrt(length(dayCT-dayST))],'-');
%             plot([8 8],[mean(nightCT-nightST)-std(nightCT-nightST)/sqrt(length(nightCT-nightST)) mean(nightCT-nightST)+std(nightCT-nightST)/sqrt(length(nightCT-nightST))],'-');
% 
figure;hold on;
subplot(121);hold on;
plot(dayST.*Dirmatrix,allST.*Dirmatrix,'ko')
subplot(122);hold on;
plot(nightST.*Dirmatrix,allST.*Dirmatrix,'ro')


%%% All data
% Stim trials
            clear LearnST;
            for i=1:6
                % 1=all, 2=day, 3=night
                    LearnST{i}{1}=(IntExp(i).deltadayFBlearning(1:end-1)+IntExp(i).deltanightFBlearning(1:end));
                    LearnST{i}{2}=IntExp(i).deltadayFBlearning(1:end);
                    LearnST{i}{3}=IntExp(i).deltanightFBlearning(1:end-1);    
            end

            allST=([LearnST{2}{1} LearnST{3}{1} LearnST{4}{1} LearnST{5}{1}]);% LearnST{6}{1}]);
            dayST=([LearnST{2}{2} LearnST{3}{2} LearnST{4}{2} LearnST{5}{2}]);% LearnST{6}{2}]);
            nightST=([LearnST{2}{3} LearnST{3}{3} LearnST{4}{3} LearnST{5}{3}]);% LearnST{6}{3}]);
% Catch trials
            clear LearnCT;
            for i=1:6
                % 1=all, 2=day, 3=night
                    LearnCT{i}{1}=(IntExp(i).deltadayCTlearning(1:end-1)+IntExp(i).deltanightCTlearning(1:end));
                    LearnCT{i}{2}=IntExp(i).deltadayCTlearning(1:end);
                    LearnCT{i}{3}=IntExp(i).deltanightCTlearning(1:end-1);    
            end

            allCT=([LearnCT{2}{1} LearnCT{3}{1} LearnCT{4}{1} LearnCT{5}{1}]);% LearnCT{6}{1}]);
            dayCT=([LearnCT{2}{2} LearnCT{3}{2} LearnCT{4}{2} LearnCT{5}{2}]);% LearnCT{6}{2}]);
            nightCT=([LearnCT{2}{3} LearnCT{3}{3} LearnCT{4}{3} LearnCT{5}{3}]);% LearnCT{6}{3}]);

figure;hold on;
subplot(121);hold on;
plot(dayST,allST,'ko')
subplot(122);hold on;
plot(nightST,allST,'ro')


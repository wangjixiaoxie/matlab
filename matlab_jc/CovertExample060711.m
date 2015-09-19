% EXAMPLES WITH APV
clearvars
            load C:\cardinal8\CovertExample_060411.mat
            % Experiment 12 - bk75bk62 - 9.27.09
                % note that APVwn are catch trials - everything else is not
                % Threshold set at 


            figure;hold on;
            mrksize=10;
                subplot(222);hold on;
                coef=-1;
                m1=round(median(Exp12.TargetingWN))-64:round(median(Exp12.TargetingWN));
                mnbase=mean(mean(Exp12.pitchACpre(m1,:)));
                %m1=260:290;
                ravgwin=20;
                % All acsf trials from that morning
                plot(timing4(fvalsACpre(1:3:end)),coef*(mean(pitchACpre(m1,1:3:end))-mnbase),'.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsACpre(1:3:end)),ravgwin),runningaverage(coef*(mean(pitchACpre(m1,1:3:end))-mnbase),ravgwin),'k','Linewidth',3)
                % Excluding first 30 minutes of apv - all notes sung
                plot(timing4(fvalsAPV(1:3:end)),coef*(mean(pitchAPV(m1,1:3:end))-mnbase),'r.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsAPV),ravgwin),runningaverage(coef*(mean(pitchAPV(m1,:))-mnbase),ravgwin),'k','Linewidth',3)
                % All notes sung with apv + wn
                plot(timing4(fvalsAPVwnCatch),coef*(mean(pitchAPVwnCatch(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsAPVwnCatch),ravgwin),runningaverage(coef*(mean(pitchAPVwnCatch(m1,:))-mnbase),ravgwin),'k','Linewidth',3)
                % All acsf trials from the next morning
                plot(timing4(fvalsACpost(34:3:146)),coef*(mean(pitchACpost(m1,34:3:146))-mnbase),'.','Markersize',mrksize)
                plot(runningaverage(timing4(fvalsACpost(34:146)),ravgwin),runningaverage(coef*(mean(pitchACpost(m1,34:146))-mnbase),ravgwin),'k','Linewidth',3)
                tpvls=timing4([fvalsACpre fvalsAPV fvalsAPVwnCatch fvalsACpost(34:146)]);
                ppvls=coef*([(mean(pitchACpre(m1,:))) (mean(pitchAPV(m1,:))) (mean(pitchAPVwnCatch(m1,:))) (mean(pitchACpost(m1,34:146)))]-mnbase);
                ravgwin=100;
                plot(runningaverage(tpvls,ravgwin),runningaverage(ppvls,ravgwin),'k','Linewidth',3)          
                plot([6505 6515],[coef*median(mean(Exp12.pitchACpost(m1,34:146))-mnbase) coef*median(mean(Exp12.pitchACpost(m1,34:146))-mnbase)],'k','Linewidth',3)    
                plot([6496 6500],[thrsh-mnbase thrsh-mnbase])
                plot([6485 6515],[0 0])
                ylim([-200 200])
                xlim([6485 6515])

                time_offset=min(timing4(fvalsAPV))-min(Exp12.timeAPVCTL);
                subplot(224);
                hold on;
                coef=-1;
                m1=250:300;
                pitchACpreCTL=pitchACpreCTL_1(:,[1:4:660]);
                timeACpreCTL=timing4(fvalsACpreCTL_1([1:4:660]));
                mnbase=mean(mean(pitchACpreCTL(m1,:)));
                %m1=round(median(Exp12.TargetingWN))-64:round(median(Exp12.TargetingWN));
                %m1=250:350;
                ravgwin=20;
                % All acsf trials from that morning
                plot(time_offset+timeACpreCTL+864.1,coef*(mean(pitchACpreCTL(m1,:))-mnbase),'.','Markersize',mrksize)
                plot(time_offset+runningaverage(timeACpreCTL+864.1,ravgwin),runningaverage(coef*(mean(pitchACpreCTL(m1,:))-mnbase),ravgwin),'k','Linewidth',3)
                % Excluding first 30 minutes of apv - all notes sung
                plot(time_offset+Exp12.timeAPVCTL,coef*(mean(Exp12.pitchAPVCTL(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(time_offset+runningaverage(Exp12.timeAPVCTL,ravgwin),runningaverage(coef*(mean(Exp12.pitchAPVCTL(m1,:))-mnbase),ravgwin),'k','Linewidth',3)
                % All notes sung with apv + wn
                plot(time_offset+Exp12.timeAPVwnCTL,coef*(mean(Exp12.pitchAPVwnCTL(m1,:))-mnbase),'r.','Markersize',mrksize)
                plot(time_offset+runningaverage(Exp12.timeAPVwnCTL,ravgwin),runningaverage(coef*(mean(Exp12.pitchAPVwnCTL(m1,:))-mnbase),ravgwin),'k','Linewidth',3)
                % All acsf trials from the next morning
                plot(time_offset+Exp12.timeACpostCTL(33:142),coef*(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase),'.','Markersize',mrksize)
                plot(time_offset+runningaverage(Exp12.timeACpostCTL(33:142),ravgwin),runningaverage(coef*(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase),ravgwin),'k','Linewidth',3)
                tpvls=time_offset+[timeACpreCTL+864.1 Exp12.timeAPVCTL Exp12.timeAPVwnCTL Exp12.timeACpostCTL(33:142)];
                ppvls=coef*([(mean(pitchACpreCTL(m1,:))) (mean(Exp12.pitchAPVCTL(m1,:))) (mean(Exp12.pitchAPVwnCTL(m1,:))) (mean(Exp12.pitchACpostCTL(m1,33:142)))]-mnbase);
                ravgwin=100;
                plot(runningaverage(tpvls,ravgwin),runningaverage(ppvls,ravgwin),'k','Linewidth',3)          
                plot([6485 6515],[0 0])
                plot([6505 6515],[coef*median(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase) coef*median(mean(Exp12.pitchACpostCTL(m1,33:142))-mnbase)],'k','Linewidth',3)
                ylim([-200 200])
                xlim([6485 6515])

   % EXAMPLES W/O APV
% Number 7
% Targeted
clearvars
load C:\cardinal8\cardinal6birds\bk76bk63\EPC7.mat
    mrksize=10;
       subplot(221);hold on;
       coef=1;
       time=150:200;%EPC7.time;
       mnbase=median(mean(pitchpreB(time,:)));
       ravgwin=30;
       plot(timing3(fvalspre_B([1:5:end])),coef*mean(pitchpreB(time,[1:5:end])-mnbase),'.','Markersize',mrksize)
       plot(timing3(fvalsWN_B([1:5:end])),coef*mean(pitchwnB(time,[1:5:end])-mnbase),'r.','Markersize',mrksize)       
       plot(timing3(fvalspostB(1:5:176)),coef*mean(pitchPostB(time,1:5:176)-mnbase),'.','Markersize',mrksize)    
       plot(runningaverage(timing3(fvalspre_B),ravgwin),coef*runningmedian(mean(pitchpreB(time,:))-mnbase,ravgwin),'k','Linewidth',3)
       plot(runningaverage(timing3(fvalsWN_B(1:5:end)),ravgwin),coef*runningmedian(mean(pitchwnB(time,1:5:end))-mnbase,ravgwin),'k','Linewidth',3)       
       plot(runningaverage(timing3(fvalspostB(1:176)),ravgwin),coef*runningmedian(mean(pitchPostB(time,1:176))-mnbase,ravgwin),'k','Linewidth',3)
       tpvls=timing3([fvalspre_B fvalsWN_B fvalspostB(1:176)]);
       ppvls=coef*([(mean(pitchpreB(time,:))) (mean(pitchwnB(time,:))) (mean(pitchPostB(time,1:176)))]-mnbase);
       ravgwin=100;
       plot(runningaverage(tpvls,ravgwin),runningaverage(ppvls,ravgwin),'k','Linewidth',3)
       hold on;plot([6870 6884],[0 0])
        plot([6882 6884],[median(mean(pitchPostB(time,1:5:176)-mnbase)) median(mean(pitchPostB(time,1:5:176)-mnbase))],'k','Linewidth',3);
        ylim([-200 200]) 
        xlim([6870 6884])
  % Non-targeted - 'A' 
      subplot(223);hold on;
       coef=1;
       time=250:350;
       mnbase=mean(mean(pitchpreA(time,:)));
       ravgwin=50;
       plot(timing3(fvalspreA([1:3:end])),coef*mean(pitchpreA(time,:)-mnbase),'.','Markersize',mrksize)
       plot(timing3(fvalswnA),coef*mean(pitchwnA(time,:)-mnbase),'r.','Markersize',mrksize)       
       plot(timing3(fvalspostA([1:3:end])),coef*mean(pitchpostA(time,:)-mnbase),'.','Markersize',mrksize)    
       plot(runningaverage(timing3(fvalspreA([1:3:end])),ravgwin),coef*runningmedian(mean(pitchpreA(time,:))-mnbase,ravgwin),'k','Linewidth',3)
       plot(runningaverage(timing3(fvalswnA),ravgwin),coef*runningmedian(mean(pitchwnA(time,:))-mnbase,ravgwin),'k','Linewidth',3)       
       plot(runningaverage(timing3(fvalspostA([1:3:end])),ravgwin),coef*runningmedian(mean(pitchpostA(time,:))-mnbase,ravgwin),'k','Linewidth',3)       
        tpvls=timing3([fvalspreA([1:3:end]) fvalswnA fvalspostA([1:3:end])]);
       ppvls=coef*([(mean(pitchpreA(time,:))) (mean(pitchwnA(time,:))) (mean(pitchpostA(time,:)))]-mnbase);
       ravgwin=100;
       plot(runningaverage(tpvls,ravgwin),runningaverage(ppvls,ravgwin),'k','Linewidth',3)
   
       hold on;plot([6870 6884],[0 0])
        plot([6882 6884],[median(mean(pitchpostA(time,:)-mnbase)) median(mean(pitchpostA(time,:)-mnbase))],'k','Linewidth',3);
        ylim([-200 200]) 
        xlim([6870 6884])
%%%

% 
% 
% 
% 
% 
%    load C:\cardinal8\CovertNTExample_060411.mat
%    close all
%    for i=7%1:length(ExperimentPC)
%        figure;hold on;
%        coef=1-2*isequal(ExperimentPC(i).DIR,'down');
%        if length(ExperimentPC(i).time)==1
%            time=ExperimentPC(i).time;
%        else
%            time=round(mean(ExperimentPC(i).time));
%        end
%        mnbase=mean(mean(ExperimentPC(i).pitchPre(time,:)));
%        plot(runningaverage(ExperimentPC(i).timePre,5),coef*runningmedian(ExperimentPC(i).pitchPre(time,:)-mnbase,5),'.')
%        plot(runningaverage(ExperimentPC(i).timeWN,5),coef*runningmedian(ExperimentPC(i).pitchWN(time,:)-mnbase,5),'r.')       
%        plot(runningaverage(ExperimentPC(i).timePost,5),coef*runningmedian(ExperimentPC(i).pitchPost(time,:)-mnbase,5),'.')
% 
%    end
%    % numbers 2,3,6,7 most promising
%    
%    % Number 2 is nice, but non-targeted notes also shift!
%                load /cardinal8/APVinRA/cardinal6birds/bk76bk63/EPC2_60511.mat
%                % EXPERIMENT 2
%                % 11.30.09 - targeted note is C - first flat stack - auto-labeled 'x'
%                        i=2;
%                        pitchWN=EPC2.pitchWN(:,[1:4:504]);
%                        timeWN=(EPC2.timeWN([1:4:504]));
% 
%                        figure;hold on;
%                        coef=1-2*isequal(EPC2.DIR,'down');
%                        if length(EPC2.time)==1
%                            time=EPC2.time;
%                        else
%                            time=round(mean(EPC2.time));
%                        end
%                        mnbase=mean(mean(EPC2.pitchPre(time,:)));
%                        ravgwin=50;
%                        plot((EPC2.timePre),coef*(EPC2.pitchPre(time,:)-mnbase),'.')
%                        plot((timeWN),coef*(pitchWN(time,:)-mnbase),'r.')
%                        plot((EPC2.timePost(1:82)),coef*(EPC2.pitchPost(time,1:82)-mnbase),'.')
% 
%                        plot(runningaverage(EPC2.timePre,ravgwin),coef*runningmedian(EPC2.pitchPre(time,:)-mnbase,ravgwin),'k','Linewidth',3)
%                        plot(runningaverage(timeWN,ravgwin),coef*runningmedian(pitchWN(time,:)-mnbase,ravgwin),'k','Linewidth',3)
%                        %plot(runningaverage(EPC2.timePost(1:82),ravgwin),coef*runningmedian(EPC2.pitchPost(time,1:82)-mnbase,ravgwin),'k','Linewidth',3)
%                        plot([8009 8012],[mean(EPC2.pitchPost(time,1:82)-mnbase) mean(EPC2.pitchPost(time,1:82)-mnbase)])
%                        plot([8002 8012],[0 0])
%                        ylim([-200 200])
%                        xlim([8002 8012])
%                 % CTL note (nt = non-targeted) is the note after 'a' - second high
%                 % stack
%                        figure;hold on;
%                        coef=1-2*isequal(EPC2.DIR,'down');
%                        time=[250:350];
%                        mnbase=mean(mean(pitchPREnt(time,:)));
%                        ravgwin=30;
%                        plot(timing4(fvalsPREnt),coef*(mean(pitchPREnt(time,:))-mnbase),'.')
%                        plot(timing4(fvalsWNnt),coef*(mean(pitchWNnt(time,:))-mnbase),'r.')
%                        plot(timing4(fvalsPOSTnt),coef*(mean(pitchPOSTnt(time,:))-mnbase),'.')
% 
%                        plot(runningaverage(timing4(fvalsPREnt),ravgwin),coef*runningmedian(mean(pitchPREnt(time,:))-mnbase,ravgwin),'k','Linewidth',3)
%                        plot(runningaverage(timing4(fvalsWNnt),ravgwin),coef*runningmedian(mean(pitchWNnt(time,:))-mnbase,ravgwin),'k','Linewidth',3)
%                        plot(runningaverage(timing4(fvalsPOSTnt),ravgwin),coef*runningmedian(mean(pitchPOSTnt(time,:))-mnbase,ravgwin),'k','Linewidth',3)
%                        plot([8032 8035],[median(mean(pitchPOSTnt(time,:))-mnbase) median(mean(pitchPOSTnt(time,:))-mnbase)])
%                        plot([8025 8035],[0 0])
%                        ylim([-200 200])
%                        xlim([8025 8035])



% 06.07.11 - Krupa


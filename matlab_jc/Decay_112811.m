load Decay2.mat
%%
% We can make a pretty solid case that learning is revealed 
%   1. within 2hrs of turning APV off
%   2. within 1hr of the first song post-APV.
% We need to compare this to the rate of positive control learning.



%% Number of songs performed
                % Positive control
                    figure; hold on;
                    avgPC=zeros(300,length(ExperimentPC));
                    count=zeros(1,300);
                    for i=1:length(ExperimentPC)
                        direction=1-2*(isequal(ExperimentPC(i).DIR,'down'));
                        for j=1:300
                            if j<length(ExperimentPC(i).timeWN)+1
                                count(j)=count(j)+1;
                                avgPC(j,count(j))=mean(direction*((ExperimentPC(i).pitchWN(ExperimentPC(i).time,j))'-mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)')))/mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)));
                            end
                        end
                    end
                    clear mnPCsyll sdPCsyll ssizePCsyll
                    for i=1:size(avgPC,1)
                        mnPCsyll(i)=mean(avgPC(i,find(avgPC(i,:)~=0)));
                        sdPCsyll(i)=std(avgPC(i,find(avgPC(i,:)~=0)));
                        ssizePCsyll(i)=length(avgPC(i,find(avgPC(i,:)~=0)));
                    end
                % APV in RA
                    figure; hold on;
                    avgAPV=zeros(300,length(ind));
                    count=zeros(1,300);
                    for z=1:length(ind)
                        i=ind(z);
                        direction=1-2*(isequal(Esave(i).DIR,'down'));
                        for j=1:300
                            if j<length(Esave(i).Tpost)+1
                                count(j)=count(j)+1;
                                avgAPV(j,count(j))=mean(direction*(Esave(i).FFpost(j)-mean(Esave(i).FFpre)))/mean(Esave(i).FFpre);
                            end
                        end
                    end
                    clear mnAPVsyll sdAPVsyll ssizeAPVsyll
                    for i=1:size(avgAPV,1)
                        mnAPVsyll(i)=mean(avgAPV(i,find(avgAPV(i,:)~=0)));
                        sdAPVsyll(i)=std(avgAPV(i,find(avgAPV(i,:)~=0)));
                        ssizeAPVsyll(i)=length(avgAPV(i,find(avgAPV(i,:)~=0)));
                    end

            figure;hold on;
            ravgwin=19;
            plot(runningaverage(mnAPVsyll,ravgwin),'r');plot(runningaverage(mnAPVsyll-sdAPVsyll./sqrt(ssizeAPVsyll*ravgwin),ravgwin),'r');plot(runningaverage(mnAPVsyll+sdAPVsyll./sqrt(ssizeAPVsyll*ravgwin),ravgwin),'r')
            plot(runningaverage(mnPCsyll,ravgwin));plot(runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot(runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
            xlim([1 100-ravgwin]) % keeps 20 of 21 APV exps and 10 of 14 CTL exps
  %%% ********************          
            figure;hold on;
            % compare first 20 songs in APV experiments with ravg=20 for CTLs
            ravgwin=19;
            g=mean(mnAPVsyll(1:20));
            s=std(mnAPVsyll(1:20))/sqrt(21);
            plot([0 140],[g g],'r');plot([0 140],[g-s g-s],'r');plot([0 140],[g+s g+s],'r')
            plot(runningaverage(mnPCsyll,ravgwin));plot(runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot(runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
            xlim([1 140]) % keeps 20 of 21 APV exps and 10 of 14 CTL exps
%%% **************************
            %%%% Are they sampled at similar rates???
            clear m n
                for i=ind
                    z=(diff(Esave(i).Tpost));
                    m(i)=median(z(find(z>0)));
                end
                for i=1:14
                    z=(diff(ExperimentPC(i).timeWN));
                    n(i)=median(z(find(z>0)));
                end
                mean(m(ind)) % 0.1545
                mean(n) % 0.1339
            % Estimate sampling rates
            clear facM facN
                for i=1:length(m)
                    if m(i)<0.1
                        facM(i)=1;
                    else
                        facM(i)=m(i)/0.1;
                    end
                end
                for i=1:length(n)
                    if  n(i)<0.1
                        facN(i)=1;
                    else
                        facN(i)=n(i)/0.1;
                    end
                end
%%%%%%%%%%%%%%%%%%

% FIG1D
ravgwin=29;
clear raAPV raPC
for i=1:size(avgAPV,2)
    raAPV(i,:)=runningaverage(avgAPV(:,i),ravgwin);
end
for i=1:size(avgPC,2)
    raPC(i,:)=runningaverage(avgPC(:,i),ravgwin);
end
figure;hold on;subplot(121);hold on;
plot(ravgwin/2+[1:7:7*length(raAPV)],mean(raAPV),'r');plot(ravgwin/2+[1:7:7*length(raAPV)],mean(raAPV)-std(raAPV)/sqrt(21),'r');plot(ravgwin/2+[1:7:7*length(raAPV)],mean(raAPV)+std(raAPV)/sqrt(21),'r')
plot(ravgwin/2+[1:7:7*length(raPC)],mean(raPC),'b');plot(ravgwin/2+[1:7:7*length(raPC)],mean(raPC)-std(raPC)/sqrt(14),'b');plot(ravgwin/2+[1:7:7*length(raPC)],mean(raPC)+std(raPC)/sqrt(14),'b')
xlim([0 300])

subplot(122);hold on;
plot([1:7:140],mean(avgAPV(1:20,:)'),'r');plot([1:7:140],mean(avgAPV(1:20,:)')+std(avgAPV(1:20,:)')/sqrt(21),'r');plot([1:7:140],mean(avgAPV(1:20,:)')-std(avgAPV(1:20,:)')/sqrt(21),'r')
plot([0 140],[0 0])


%%%%%%%%%%%%%%%%%%
%%% FIG1A             
figure;hold on;subplot(121);hold on;
ravgwin=24;
plot(ravgwin/2+[1:7:length(runningaverage(mnAPVsyll,ravgwin))*7],100*runningaverage(mnAPVsyll,ravgwin),'r');plot(ravgwin/2+[1:7:length(runningaverage(mnAPVsyll,ravgwin))*7],100*runningaverage(mnAPVsyll-sdAPVsyll./sqrt(ssizeAPVsyll*ravgwin),ravgwin),'r');plot(ravgwin/2+[1:7:length(runningaverage(mnAPVsyll,ravgwin))*7],100*runningaverage(mnAPVsyll+sdAPVsyll./sqrt(ssizeAPVsyll*ravgwin),ravgwin),'r')
plot(ravgwin/2+[1:7:length(runningaverage(mnPCsyll,ravgwin))*7],100*runningaverage(mnPCsyll,ravgwin));plot(ravgwin/2+[1:7:length(runningaverage(mnPCsyll,ravgwin))*7],100*runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot(ravgwin/2+[1:7:length(runningaverage(mnPCsyll,ravgwin))*7],100*runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
xlim([1 300]) % keeps 20 of 21 APV exps and 10 of 14 CTL exps
ylim([0 1.2])
subplot(122);hold on;
ravgwin=29;
plot(ravgwin/2+[1:7:length(runningaverage(mnAPVsyll,ravgwin))*7],100*runningaverage(mnAPVsyll,ravgwin),'r');plot(ravgwin/2+[1:7:length(runningaverage(mnAPVsyll,ravgwin))*7],100*runningaverage(mnAPVsyll-sdAPVsyll./sqrt(ssizeAPVsyll*ravgwin),ravgwin),'r');plot(ravgwin/2+[1:7:length(runningaverage(mnAPVsyll,ravgwin))*7],100*runningaverage(mnAPVsyll+sdAPVsyll./sqrt(ssizeAPVsyll*ravgwin),ravgwin),'r')
plot(ravgwin/2+[1:7:length(runningaverage(mnPCsyll,ravgwin))*7],100*runningaverage(mnPCsyll,ravgwin));plot(ravgwin/2+[1:7:length(runningaverage(mnPCsyll,ravgwin))*7],100*runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot(ravgwin/2+[1:7:length(runningaverage(mnPCsyll,ravgwin))*7],100*runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
xlim([1 300]) % keeps 20 of 21 APV exps and 10 of 14 CTL exps
ylim([0 1.2])
%
figure;hold on;
plot(mnAPVsyll(1:10),'r');plot(mnAPVsyll(1:10)+sdAPVsyll(1:10)./sqrt(ssizeAPVsyll(1:10)),'r');plot(mnAPVsyll(1:10)-sdAPVsyll(1:10)./sqrt(ssizeAPVsyll(1:10)),'r');
plot(mnPCsyll(1:10),'b');plot(mnPCsyll(1:10)+sdPCsyll(1:10)./sqrt(ssizePCsyll(1:10)),'b');plot(mnPCsyll(1:10)-sdPCsyll(1:10)./sqrt(ssizePCsyll(1:10)),'b');

xlim([1 300]) % keeps 20 of 21 APV exps and 10 of 14 CTL exps
ylim([0 1.2])

  %%% ********************          
  %%% FIG1B
            figure;hold on;
            % compare first 20 songs in APV experiments with ravg=20 for CTLs
            ravgwin=24;
            g=mean(mnAPVsyll(1:25));
            s=std(mnAPVsyll(1:25))/sqrt(21);
            plot([0 700],[g g],'r');plot([0 700],[g-s g-s],'r');plot([0 700],[g+s g+s],'r')
            plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll,ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
            xlim([1 600]) % keeps 21 of 21 APV exps and 10 of 14 CTL exps
      
%%% **************************
%%% FIG1C
    % Adjust this to only include experiments with CV recovery on first
    % song post-apv
    indX=[3:5 8 11:12 14:21]; % n=14
                    clear mnAPVsyll2 sdAPVsyll2 ssizeAPVsyll2
                    avgAPV2=avgAPV(:,indX);
                    for i=1:size(avgAPV2,1)
                        mnAPVsyll2(i)=mean(avgAPV2(i,find(avgAPV2(i,:)~=0)));
                        sdAPVsyll2(i)=std(avgAPV2(i,find(avgAPV2(i,:)~=0)));
                        ssizeAPVsyll2(i)=length(avgAPV2(i,find(avgAPV2(i,:)~=0)));
                    end
            figure;hold on;
            % compare first 20 songs in APV experiments with ravg=20 for CTLs
            ravgwin=24;
            g=mean(mnAPVsyll2(1:25));
            s=std(mnAPVsyll2(1:25))/sqrt(length(indX));
            plot([0 700],[g g],'r');plot([0 700],[g-s g-s],'r');plot([0 700],[g+s g+s],'r')
            plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll,ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
            xlim([1 600]) % keeps 21 of 21 APV exps and 10 of 14 CTL exps

%%%%%




%% Amount of time singing --- interpolation
                %APV
                    clear INTvals
                    for z=1:length(ind)
                        i=ind(z);
                        direction=1-2*(isequal(Esave(i).DIR,'down'));
                        vls=find(Esave(i).Tpost-min(Esave(i).Tpost)<5);
                        [a,b]=unique(Esave(i).Tpost(vls));
                        b=[1 b];
                        X=[];Y=[];
                        for j=2:length(b)
                            X=[X Esave(i).Tpost(vls(b(j-1)+1):vls(b(j)))-min(Esave(i).Tpost)+0.001*[1:1:b(j)-b(j-1)]];
                            Y=[Y direction*((Esave(i).FFpost(vls(b(j-1)+1):vls(b(j))))-mean(Esave(i).FFpre))/mean(Esave(i).FFpre)];
                        end
                        INTvals(i,:)=interp1(X,Y,[0:0.05:5],'linear');
                    end
                    clear mnINT sdINT ssizeINT
                    INTvals=INTvals(ind,:);
                    for k=1:size(INTvals,2)
                            mnINT(k)=mean(INTvals(~isnan(INTvals(:,k)),k));
                            sdINT(k)=std(INTvals(~isnan(INTvals(:,k)),k));
                            ssizeINT(k)=length(INTvals(~isnan(INTvals(:,k)),k));
                    end
                % CTL
                clear EsaveCTL
                for i=1:14
                    if length(ExperimentPC(i).time)>1
                        EsaveCTL(i).FFwn=mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,:));
                        EsaveCTL(i).FFpre=mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:));
                    else
                        EsaveCTL(i).FFwn=(ExperimentPC(i).pitchWN(ExperimentPC(i).time,:));
                        EsaveCTL(i).FFpre=(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:));
                    end
                    EsaveCTL(i).Tpre=ExperimentPC(i).timePre;
                    EsaveCTL(i).Twn=ExperimentPC(i).timeWN;
                    EsaveCTL(i).DIR=ExperimentPC(i).DIR;
                end
                %
                clear INTvalsCTL
                for i=1:14
                    direction=1-2*(isequal(EsaveCTL(i).DIR,'down'));
                    vls=find(EsaveCTL(i).Twn-min(EsaveCTL(i).Twn)<5);
                    [a,b]=unique(EsaveCTL(i).Twn(vls));
                    b=[1 b];
                    X=[];Y=[];
                    for j=2:length(b)
                        X=[X EsaveCTL(i).Twn(vls(b(j-1)+1):vls(b(j)))-min(EsaveCTL(i).Twn)+0.001*[1:1:b(j)-b(j-1)]];
                        Y=[Y direction*((EsaveCTL(i).FFwn(vls(b(j-1)+1):vls(b(j))))-mean(EsaveCTL(i).FFpre))/mean(EsaveCTL(i).FFpre)];
                    end
                    INTvalsCTL(i,:)=interp1(X,Y,[0:0.05:5],'linear');
                end
                clear mnINTCTL sdINTCTL ssizeINTCTL
                for k=1:size(INTvalsCTL,2)
                        mnINTCTL(k)=mean(INTvalsCTL(~isnan(INTvalsCTL(:,k)),k));
                        sdINTCTL(k)=std(INTvalsCTL(~isnan(INTvalsCTL(:,k)),k));
                        ssizeINTCTL(k)=length(INTvalsCTL(~isnan(INTvalsCTL(:,k)),k));
                end
%%% FIG2A
    figure;hold on;
    plot(mnINT,'r');plot(mnINT-sdINT./sqrt(ssizeINT),'r');plot(mnINT+sdINT./sqrt(ssizeINT),'r')
    plot(mnINTCTL);plot(mnINTCTL-sdINTCTL./sqrt(ssizeINTCTL));plot(mnINTCTL+sdINTCTL./sqrt(ssizeINT))
%%% FIG2B 
    figure;hold on;
    ravgwin=9;
    plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINT,ravgwin),'r','Linewidth',3);plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINT-sdINT./sqrt(ssizeINT),ravgwin),'r');plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINT+sdINT./sqrt(ssizeINT),ravgwin),'r')
    plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINTCTL,ravgwin),'Linewidth',3);plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINTCTL-sdINTCTL./sqrt(ssizeINTCTL),ravgwin));plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINTCTL+sdINTCTL./sqrt(ssizeINT),ravgwin))
    xlim([0 3.3]) % sample size drops off dramatically after this time
    plot([0 3.3],[0 0])
    ylim([-0.01 0.018])
%%% FIG2C
        figure;hold on;
        % 0.0083 and 0.0029 were calculated from the code below ---> They are the
            % first half hour of singing with ACSF in RA 
        plot([0 3.5],[0.0083 0.0083],'r');plot([0 3.5],[0.0083-0.0029 0.0083-0.0029],'r');plot([0 3.5],[0.0083+0.0029 0.0083+0.0029],'r')
        ravgwin=4;
        plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINTCTL,ravgwin));plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINTCTL-sdINTCTL./sqrt(ssizeINTCTL),ravgwin));plot(runningaverage([0:0.05:5],ravgwin),runningaverage(mnINTCTL+sdINTCTL./sqrt(ssizeINT),ravgwin))
        xlim([0 3.5]) % sample size drops off dramatically after this time
        plot([0 3.5],[0 0])



%% Earlier analysis - group per hour (could be adjusted for 30min, 2hrs, etc)


        % Hours after APVwn off 
                clear mnhr
                figure;hold on;
                for z=1:length(ind)
                    i=ind(z);
                    direction=1-2*(isequal(Esave(i).DIR,'down'));
                    plot(Esave(i).Tpost-Esave(i).Tapvend,direction*(Esave(i).FFpost-mean(Esave(i).FFpre)),'.')
                    for k=1:10
                        ishr=find(k-1<Esave(i).Tpost-Esave(i).Tapvend & Esave(i).Tpost-Esave(i).Tapvend<k);
                        mnhr(z,k)=mean(direction*(Esave(i).FFpost(ishr)-mean(Esave(i).FFpre)))/mean(Esave(i).FFpre);
                    end
                end
                mnhrA=mnhr;
                clear mnhr2 ssize sdhr
                figure;plot(mnhr')
                for i=1:10
                    mnhr2A(i)=mean(mnhr(find(~isnan(mnhr(:,i))),i));
                    sdhrA(i)=std(mnhr(find(~isnan(mnhr(:,i))),i));
                    ssizeA(i)=length(find(~isnan(mnhr(:,i))));
                end
                figure;hold on;
                errorbar(mnhr2A,sdhrA./sqrt(ssizeA))
                plot(mnhr','.')

        % Hours after first song Post WN
                clear mnhr
                for z=1:length(ind)
                    i=ind(z);
                    direction=1-2*(isequal(Esave(i).DIR,'down'));
                    %plot(Esave(i).Tpost-min(Esave(i).Tpost),direction*(Esave(i).FFpost-mean(Esave(i).FFpre)),'.')
                    for k=1:10
                        ishr=find(k-1<Esave(i).Tpost-min(Esave(i).Tpost) & Esave(i).Tpost-min(Esave(i).Tpost)<k);
                        mnhr(z,k)=mean(direction*(Esave(i).FFpost(ishr)-mean(Esave(i).FFpre)))/mean(Esave(i).FFpre);
                    end
                end
                mnhrB=mnhr;
                clear mnhr2 ssize sdhr
                figure;plot(mnhr')
                for i=1:10
                    mnhr2B(i)=mean(mnhr(find(~isnan(mnhr(:,i))),i));
                    sdhrB(i)=std(mnhr(find(~isnan(mnhr(:,i))),i));
                    ssizeB(i)=length(find(~isnan(mnhr(:,i))));
                end
                figure;hold on;
                errorbar(mnhr2B,sdhrB./sqrt(ssizeB))
                plot(mnhr','.')
        % Compare to Positive control experiments (normal progression of learning)
        figure;hold on;
                    clear mnhr
                for i=1:length(ExperimentPC)
                    direction=1-2*(isequal(ExperimentPC(i).DIR,'down'));
                    for k=1:10
                        ishr=find(k-1<ExperimentPC(i).timeWN-min(ExperimentPC(i).timeWN) & ExperimentPC(i).timeWN-min(ExperimentPC(i).timeWN)<k);
                        mnhr(i,k)=mean(direction*(mean(ExperimentPC(i).pitchWN(ExperimentPC(i).time,ishr)')-mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)')))/mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)));
                    end
                end
                mnhrPC=mnhr;
                clear mnhr2PC ssizePC sdhrPC
                figure;plot(mnhr')
                for i=1:10
                    mnhr2PC(i)=mean(mnhr(find(~isnan(mnhr(:,i))),i));
                    sdhrPC(i)=std(mnhr(find(~isnan(mnhr(:,i))),i));
                    ssizePC(i)=length(find(~isnan(mnhr(:,i))));
                end
                figure;hold on;
                errorbar(mnhr2PC,sdhrPC./sqrt(ssizePC))
                plot(mnhr','.')

                figure;hold on;
                errorbar(mnhr2PC,sdhrPC./sqrt(ssizePC))
                errorbar(mnhr2A,sdhrA./sqrt(ssizeA),'r')
                errorbar(mnhr2B,sdhrB./sqrt(ssizeB),'g')

        
%% 11.29.11 4:20pm - helping FIG1C
% 'slices'
load Workspace112911.mat
%%% FIG1C 
    %  experiments with total CV recovery on first song post-apv
    indX=[3:5 8 11:12 14:21]; % n=14
                    clear mnAPVsyll2 sdAPVsyll2 ssizeAPVsyll2
                    avgAPV2=avgAPV(:,indX);
                    for i=1:size(avgAPV2,1)
                        mnAPVsyll2(i)=mean(avgAPV2(i,find(avgAPV2(i,:)~=0)));
                        sdAPVsyll2(i)=std(avgAPV2(i,find(avgAPV2(i,:)~=0)));
                        ssizeAPVsyll2(i)=length(avgAPV2(i,find(avgAPV2(i,:)~=0)));
                    end
            figure;hold on;
            % compare first 25 songs in APV experiments with ravg=25 for CTLs
            ravgwin=24;
            g=mean(mnAPVsyll2(1:25));
            s=std(mnAPVsyll2(1:25))/sqrt(length(indX));
            plot([0 700],[g g],'r');plot([0 700],[g-s g-s],'r');plot([0 700],[g+s g+s],'r')
            plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll,ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
            xlim([1 600]) % keeps 21 of 21 APV exps and 10 of 14 CTL exps
%%% FIG1D
%  experiments with a small amount of time (<3hrs) before the first song
% performance post-APV (note that intersection between this and indX is
% only n=3
indY=[1 2 5 6 7 8 10 12 13 16];         
clear mnAPVsyll3 sdAPVsyll3 ssizeAPVsyll3
                    avgAPV3=avgAPV(:,indY);
                    for i=1:size(avgAPV3,1)
                        mnAPVsyll3(i)=mean(avgAPV3(i,find(avgAPV3(i,:)~=0)));
                        sdAPVsyll3(i)=std(avgAPV3(i,find(avgAPV3(i,:)~=0)));
                        ssizeAPVsyll3(i)=length(avgAPV3(i,find(avgAPV3(i,:)~=0)));
                    end
            figure;hold on;
            % compare first 30 songs in APV experiments with ravg=30 for CTLs
            ravgwin=24;
            g=mean(mnAPVsyll3(1:25));
            s=std(mnAPVsyll3(1:25))/sqrt(length(indY));
            plot([0 700],[g g],'r');plot([0 700],[g-s g-s],'r');plot([0 700],[g+s g+s],'r')
            plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll,ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll-sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin));plot([1:7:length(runningaverage(mnPCsyll,ravgwin))*7],runningaverage(mnPCsyll+sdPCsyll./sqrt(ssizePCsyll*ravgwin),ravgwin))
            xlim([1 600]) % keeps 21 of 21 APV exps and 10 of 14 CTL exps




% These are the experiments where we allowed the birds to sing shortly after having APV in RA
    for i=ind%[1 2 5 6 7 8 10 12 13 16]
        SDpre=std(Esave(i).FFpre(end-50:end));
        figure;hold on;
        plot([0 100],[SDpre SDpre],'b')
        plot(runningstd(Esave(i).FFpost,20),'r')
    end
% Among these, 5, 6, 7, 8, 10, 12 show fairly rapid CV recovery



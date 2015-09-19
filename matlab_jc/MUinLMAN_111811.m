load /bulbul5/TWmusc/DSumm6.mat
% load DSumm4.mat
% Note that all .deltaFF... fields are normalized for learning direction
% e.g. .deltaFF already has coef taken into account
% e.g. .PosCTL_deltaFFperwakinghour also has coef taken into account

% 1 - pk20r49, Exp1
% 2 - pk20r49, Exp2
% 3 - pk20r49, Exp3

% 4 - pu34, Exp1 
% 5 - pu34, Exp2

% 6 - bk20bk45, Exp1, nt 'a' (upshift)
% 7 - bk20bk45, Exp1, nt 'b' (downshift - how much did he learn here?)
% 8 - bk20bk45, Exp2, nt 'a' 
% 9 - bk20bk45, Exp2, nt 'b' 
% 10 - bk20bk45, Exp3, nt 'a' 
% 11 - bk20bk45, Exp3, nt 'b' 
% 12 - bk20bk45, Exp5, nt 'a' (Exp4 DUR had no escapes for either 'a' or 'b')
% 13 - bk20bk45, Exp5, nt 'b' 

% 14 - bk28w6, Exp1

%% What about during MU infusion?
clear delFFmu delFFmu2
figure;hold on;
for i=1:14
    subplot(4,4,i);hold on;
    Tkeep=find(timing3(DSumm(i).fv{2})>min(timing3(DSumm(i).fv{2}))+0.5);
    FF=DSumm(i).coef*(mean(DSumm(i).pitch{2}(DSumm(i).window,Tkeep))'-DSumm(i).mnFFpre);
    indkeep=find(FF>median(FF)-1.5*(prctile(FF,75)-prctile(FF,25)) & FF<median(FF)+1.5*(prctile(FF,75)-prctile(FF,25)));
    FF=FF(indkeep);
    FFcv=mean(DSumm(i).pitch{2}(DSumm(i).window,Tkeep));
    FFcv=FFcv(indkeep);
    CV(i)=std(FFcv)/mean(FFcv);
    plot(FF)
    X=[ones(length(FF),1) [1:1:length(FF)]'];
    Y=FF;
    [b,bint] = regress(Y,X);
    delFFmu(i)=b(2)*size(Y,1)/DSumm(i).mnFFpre;
    %delFFmu2(i)=DSumm(i).coef*(median(FF(1:10))-median(FF(end-10:end)));
end

mean(delFFmu) % -0.19
std(delFFmu)/sqrt(14) % 0.37

for i=1:14;CVpre(i)=DSumm(i).CVpre;end
mean(CV)/mean(CVpre)
std(CV./CVpre)/sqrt(14)
%%

%%%% Look at examples --- 12.04.11 %%%%

% exclude first 30min
DSumm(9).fv{2}=DSumm(9).fv{2}(26:end);
DSumm(9).pitch{2}=DSumm(9).pitch{2}(:,26:end);
DSumm(9).istrig{2}=DSumm(9).istrig{2}(26:end);

figure;hold on;
num=[4 9 13];
num=4
ravgwin=25;
for z=1:length(num)
    i=num(z);
    subplot(1,3,z);hold on;
    for j=1:3
        if 1
            plot(timing3(DSumm(i).fv{j}(find(DSumm(i).istrig{j})))-min(timing3(DSumm(i).fv{2})),(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,(find(DSumm(i).istrig{j}))))-DSumm(i).mnFFpre)),'r.')
            plot(timing3(DSumm(i).fv{j}(find(~DSumm(i).istrig{j})))-min(timing3(DSumm(i).fv{2})),(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,(find(~DSumm(i).istrig{j}))))-DSumm(i).mnFFpre)),'k.')
           
        else
            plot(timing3(DSumm(i).fv{j})-min(timing3(DSumm(i).fv{2})),(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre)),'b.')
        end
       plot([min(timing3(DSumm(i).fv{j})-min(timing3(DSumm(i).fv{2}))) max(timing3(DSumm(i).fv{j})-min(timing3(DSumm(i).fv{2})))],[mean(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre)) mean(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre))],'Linewidth',3,'color','k')
       plot(runningaverage(timing3(DSumm(i).fv{j}),ravgwin)-min(timing3(DSumm(i).fv{2})),runningaverage(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre),ravgwin),'g')     
    end
    ylim([-150 150])
    xlim([-2 7])
end
%% Best
figure;hold on;
num=4;
ravgwin=75;
for z=1:length(num)
    i=num(z);
    %subplot(4,4,z);hold on;
    for j=1:3
        if j==2
            plot(timing3(DSumm(i).fv{j}(find(DSumm(i).istrig{j})))-min(timing3(DSumm(i).fv{2})),(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,(find(DSumm(i).istrig{j}))))-DSumm(i).mnFFpre)),'r.')
            plot(timing3(DSumm(i).fv{j}(find(~DSumm(i).istrig{j})))-min(timing3(DSumm(i).fv{2})),(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,(find(~DSumm(i).istrig{j}))))-DSumm(i).mnFFpre)),'k.')
           
        else
            plot(timing3(DSumm(i).fv{j})-min(timing3(DSumm(i).fv{2})),(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre)),'b.')
        end
       plot([min(timing3(DSumm(i).fv{j})-min(timing3(DSumm(i).fv{2}))) max(timing3(DSumm(i).fv{j})-min(timing3(DSumm(i).fv{2})))],[mean(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre)) mean(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre))],'Linewidth',3,'color','k')
       plot(runningaverage(timing3(DSumm(i).fv{j}),ravgwin)-min(timing3(DSumm(i).fv{2})),runningaverage(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre),ravgwin),'g')     
    end
    ylim([-150 150])
    xlim([-2 5.5])
end
%% Control example
load /bulbul4/CovertAnalysis/TWMUctl.mat
figure;hold on;
plot(runningaverage(24*(Tvals(1:28)-min(Tvals(1:28))),5),runningaverage(Fvals(1:28),5)-mean(Fvals(1:10)),'-')
plot([0 4.5],[0 0],'-')
xlim([0 4.5]);ylim([-60 60])


%% Best
figure;hold on;
num=4;
ravgwin=50;
clear tkeep fkeep
    i=num;
    %subplot(4,4,z);hold on;
    for j=1:3
       tkeep{j}=runningaverage(timing3(DSumm(i).fv{j}),ravgwin)-min(timing3(DSumm(i).fv{2}));
       fkeep{j}=runningmedian(DSumm(i).coef*(mean(DSumm(i).pitch{j}(DSumm(i).window,:))-DSumm(i).mnFFpre),ravgwin);
    end
    plot([tkeep{1}+1 tkeep{2} ],[fkeep{1} fkeep{2} ],'b')
    plot(tkeep{3},fkeep{3},'b')
    plot(tkeep{2},fkeep{2},'r')
    plot([-0.5 6],[0 0])
    ylim([-60 60])
    xlim([-0.5 6])
%% Control
figure;hold on
mnpre=mean(mean(pitchPre(DSumm(4).twindow,1:20)));
plot(runningaverage(timing3(fv4(1:120)),20)-min(timing3(fv4)),runningaverage(mean(pitchPre(DSumm(4).twindow,1:120)),20)-mnpre,'r')
    plot([-0.5 4],[0 0])
    ylim([-60 60])
    xlim([-0.5 4])

    
    
%%
%%%% Lab meeting response controls --- 11.30.11 %%%%%%%%

%%%%% MATCHED PER SHIFT -normalized per hour - grouped by shift

    load /bulbul5/TWmusc/DSumm4.mat
    clear ctlshift actshift
    for i=1:14
        ctlshift(i)=(DSumm(i).PosCTL_deltaFFperwakinghour)/DSumm(i).mnFFpre;
        actshift(i)=(DSumm(i).deltaFF/DSumm(i).musctime)/DSumm(i).mnFFpre;
    end
    ctlshifts(1)=mean(ctlshift([1:3]));ctlshifts(2)=mean(ctlshift([4:5]));ctlshifts(3)=mean(ctlshift([6 8 10 12]));ctlshifts(4)=mean(ctlshift([7 9 11 13]));ctlshifts(5)=ctlshift(14);
    actshifts(1)=mean(actshift([1:3]));actshifts(2)=mean(actshift([4:5]));actshifts(3)=mean(actshift([6 8 10 12]));actshifts(4)=mean(actshift([7 9 11 13]));actshifts(5)=actshift(14);
    
        figure;plot(100*ctlshifts,100*actshifts,'.','Markersize',15);hold on;plot(100*[-0.005 0.005],100*[-0.005 0.005])
        plot(100*[mean(ctlshifts)-std(ctlshifts)/sqrt(length(ctlshifts)) mean(ctlshifts)+std(ctlshifts)/sqrt(length(ctlshifts))],100*[mean(actshifts) mean(actshifts)],'r')
        plot(100*[mean(ctlshifts) mean(ctlshifts)],100*[mean(actshifts)-std(actshifts)/sqrt(length(actshifts)) mean(actshifts)+std(actshifts)/sqrt(length(actshifts))],'r')
        xlim(100*[-0.005 0.005]);ylim(100*[-0.005 0.005])
        plot([-0.5 0.5],[0 0]);plot([0 0],[-0.5 0.5])

%%%%% MATCHED PER SHIFT -normalized per hour
    load /bulbul5/TWmusc/DSumm4.mat
    clear ctlshift actshift
    for i=1:14
        ctlshift(i)=(DSumm(i).PosCTL_deltaFFperwakinghour)/DSumm(i).mnFFpre;
        actshift(i)=(DSumm(i).deltaFF/DSumm(i).musctime)/DSumm(i).mnFFpre;
    end
        figure;plot(100*ctlshift,100*actshift,'.','Markersize',15);hold on;plot(100*[-0.005 0.005],100*[-0.005 0.005])
        plot(100*[mean(ctlshift)-std(ctlshift)/sqrt(length(ctlshift)) mean(ctlshift)+std(ctlshift)/sqrt(length(ctlshift))],100*[mean(actshift) mean(actshift)],'r')
        plot(100*[mean(ctlshift) mean(ctlshift)],100*[mean(actshift)-std(actshift)/sqrt(length(actshift)) mean(actshift)+std(actshift)/sqrt(length(actshift))],'r')
        
        xlim(100*[-0.005 0.005]);ylim(100*[-0.005 0.005])
        plot([-0.5 0.5],[0 0]);plot([0 0],[-0.5 0.5])

%%%%% MATCHED PER DAY (i.e. PER INACTIVATION)
        load /bulbul5/TWmusc/DSumm6.mat
        load /bulbul5/TWmusc/ptAM2.mat
        clear bpre bintpre ppre delTpre
        for i=1:length(ptAM)
            X=[ones(length(tAM{i}),1) tAM{i}'-min(tAM{i})];
            delTpre(i)=max(tAM{i})-min(tAM{i});
            Y=ptAM{i}'-min(ptAM{i});
            [b,bint] = regress(Y,X);
            bpre(i)=DSumm(i).coef*b(2);
            bintpre(i,:)=DSumm(i).coef*bint(2,:);
            ppre(i,:)=polyfit(X(:,2),Y,1);
        end
        clear bpost bintpost ppost delTpost
        for i=1:length(ptPM)
            X=[ones(length(tPM{i}),1) tPM{i}'-min(tPM{i})];
            delTpost(i)=max(tPM{i})-min(tPM{i});
            Y=ptPM{i}'-min(ptPM{i});
            [b,bint] = regress(Y,X);
            bpost(i)=DSumm(i).coef*b(2);
            bintpost(i,:)=DSumm(i).coef*bint(2,:);
            ppost(i,:)=polyfit(X(:,2),Y,1);
        end

        % Matched comparison
        
        clear ctlshift actshift
        for i=1:14
            muschrs(i)=DSumm(i).musctime;
            ctlshift(i)=muschrs(i)*((((delTpre(i)*bpre(i)+delTpost(i)*bpost(i))/(delTpre(i)+delTpost(i))))/DSumm(i).mnFFpre); 
            actshift(i)=muschrs(i)*((DSumm(i).deltaFF/DSumm(i).musctime)/DSumm(i).mnFFpre);
        end
        figure;plot(100*ctlshift,100*actshift,'.','Markersize',15);hold on;plot(100*[-0.005 0.015],100*[-0.005 0.015])
       % plot(100*mean(ctlshift),100*mean(actshift),'r.','Markersize',25)
        plot(100*[mean(ctlshift)-std(ctlshift)/sqrt(length(ctlshift)) mean(ctlshift)+std(ctlshift)/sqrt(length(ctlshift))],100*[mean(actshift) mean(actshift)],'r')
        plot(100*[mean(ctlshift) mean(ctlshift)],100*[mean(actshift)-std(actshift)/sqrt(length(actshift)) mean(actshift)+std(actshift)/sqrt(length(actshift))],'r')
        xlim(100*[-0.005 0.015]);ylim(100*[-0.005 0.015])
        plot([-0.5 1.5],[0 0]);plot([0 0],[-0.5 1.5])
 
        [h,p]=signtest(ctlshift-actshift)
        
        
    for i=1:length(DSumm);normfac(i)=DSumm(i).mnFFpre;end
% 1. Positive control: average amount of learning acquired during the day 
    % during that pitch shift (before asymptote) in a matched # of hours 
    
load /bulbul5/TWmusc/Data030712.mat    
figure; hold on;
    bar(1,mean(ctlshift))
    plot(1,ctlshift,'.','markersize',15)
    plot([1 1],[mean(ctlshift)-std(ctlshift)/sqrt(length(DSumm)) mean(ctlshift)+std(ctlshift)/sqrt(length(DSumm))],'Linewidth',2)

    % 2. Expression of learning with MU in LMAN
    bar(2.5,mean(delFFmu))
    plot([2.5 2.5],[mean(delFFmu)-std(delFFmu)/sqrt(length(DSumm)) mean(delFFmu)+std(delFFmu)/sqrt(length(DSumm))],'Linewidth',2)
    plot(2.5,delFFmu,'.','markersize',15)
    % 3. Acquisition of learning with MU in LMAN
    for i=1:length(DSumm)
        ACQ(i)=DSumm(i).deltaFF;
    end
    mean(ACQ./normfac) % 0.0007 (i.e. 0.07%) - nice
    std(ACQ./normfac)/sqrt(length(DSumm)) % 0.002 (i.e. 0.2%) - nice
    bar(4,mean(ACQ./normfac))
    plot([4 4],[mean(ACQ./normfac)-std(ACQ./normfac)/sqrt(length(DSumm)) mean(ACQ./normfac)+std(ACQ./normfac)/sqrt(length(DSumm))],'Linewidth',2)
    plot(4,ACQ./normfac,'.','markersize',15)
    
        
        
%%%%%
load /bulbul5/TWmusc/DSumm4.mat
        % compare slope before and after muscimol to judge if washout has occurred
        clear bpre bintpre bpost bintpost
        for i=1:length(DSumm)
            Y=mean(DSumm(i).pitch{1}(DSumm(i).window,:))-min(mean(DSumm(i).pitch{1}(DSumm(i).window,:)));
            Y=Y';
            X=[ones(length(DSumm(i).fv{1}),1) (timing3(DSumm(i).fv{1})-min(timing3(DSumm(i).fv{1})))'];
            [b,bint] = regress(Y,X);
            bpre(i)=DSumm(i).coef*b(2);
            bintpre(i,:)=DSumm(i).coef*bint(2,:);
            Y=mean(DSumm(i).pitch{3}(DSumm(i).window,:))-min(mean(DSumm(i).pitch{3}(DSumm(i).window,:)));
            Y=Y';
            X=[ones(length(DSumm(i).fv{3}),1) (timing3(DSumm(i).fv{3})-min(timing3(DSumm(i).fv{3})))'];
            [b,bint] = regress(Y,X);
            bpost(i)=DSumm(i).coef*b(2);
            bintpost(i,:)=DSumm(i).coef*bint(2,:);
        end
        % All bintpost overlap zero, indicating that none are significant
        % Turns out this isn't such a great analysis - sample size not big
        % enough
        
        
%%%% ANALYSES %%%%%%%
load /bulbul5/TWmusc/DSumm4.mat
% normalize by mean FF
    for i=1:length(DSumm);normfac(i)=DSumm(i).mnFFpre;end
% 1. Positive control: average amount of learning acquired during the day 
    % during that pitch shift (before asymptote) in a matched # of hours 
    for i=1:length(DSumm)
        pCTL(i)=DSumm(i).PosCTL_deltaFFperwakinghour*DSumm(i).musctime;
    end
    mean(pCTL./normfac) % 0.009 (i.e. 0.9%) - nice
    std(pCTL./normfac)/sqrt(length(DSumm)) %0.0009 (i.e. 0.09%) - nice
figure; hold on;
    bar(1,mean(pCTL./normfac))
    plot(1,pCTL./normfac,'.','markersize',15)
    plot([1 1],[mean(pCTL./normfac)-std(pCTL./normfac)/sqrt(length(DSumm)) mean(pCTL./normfac)+std(pCTL./normfac)/sqrt(length(DSumm))],'Linewidth',2)

    % 2. Expression of learning with MU in LMAN
    bar(2.5,mean(delFFmu))
    plot([2.5 2.5],[mean(delFFmu)-std(delFFmu)/sqrt(length(DSumm)) mean(delFFmu)+std(delFFmu)/sqrt(length(DSumm))],'Linewidth',2)
    plot(2.5,delFFmu,'.','markersize',15)
    % 3. Acquisition of learning with MU in LMAN
    for i=1:length(DSumm)
        ACQ(i)=DSumm(i).deltaFF;
    end
    mean(ACQ./normfac) % 0.0007 (i.e. 0.07%) - nice
    std(ACQ./normfac)/sqrt(length(DSumm)) % 0.002 (i.e. 0.2%) - nice
    bar(4,mean(ACQ./normfac))
    plot([4 4],[mean(ACQ./normfac)-std(ACQ./normfac)/sqrt(length(DSumm)) mean(ACQ./normfac)+std(ACQ./normfac)/sqrt(length(DSumm))],'Linewidth',2)
    plot(4,ACQ./normfac,'.','markersize',15)
    
    
 % 3. CV recovery? #13 is questionable
    for i=1:length(DSumm)
        CVrecov(i)=DSumm(i).CVpost/DSumm(i).CVpre;
    end
    figure;hold on;
    bar(5.5,1)
    bar(6.5,mean(CV./CVpre))    
    plot(6.5,CV./CVpre,'o','markersize',15)
    bar(7.5,mean(CVrecov))
    plot(7.5,CVrecov,'o','markersize',15)   
    plot([7.5 7.5],[mean(CVrecov)-std(CVrecov)/sqrt(length(CVrecov)) mean(CVrecov)+std(CVrecov)/sqrt(length(CVrecov))],'Linewidth',2)
    plot([6.5 6.5],[mean(CV./CVpre)-std(CV./CVpre)/sqrt(length(CV./CVpre)) mean(CV./CVpre)+std(CV./CVpre)/sqrt(length(CV./CVpre))],'Linewidth',2)
    xlim([0 4])
    ylim([0 1.1])
    

    figure;plot(CVrecov,'.');hold on;plot([0 15],[1 1])
% 3. Hit rate? #4 is questionable (note that I already excluded bk20bk45 exp4 a & b)
    for i=1:length(DSumm)
        Hitrate(i)=DSumm(i).HitRate;
    end
    figure;plot(Hitrate,'.');ylim([0 1]);xlim([0 15])
% 4. Expression of learning with MU in LMAN
% *********  I would need to label Exps 1-9 DUR non-catch trials to calculate this
% excluding #4 (permanently) and #13 (for now)
    figure; hold on;
    ind2=[1:3 5:12 14];
        bar(1,mean(pCTL(ind2)./normfac(ind2)))
        plot([1 1],[mean(pCTL(ind2)./normfac(ind2))-std(pCTL(ind2)./normfac(ind2))/sqrt(length(DSumm)-2) mean(pCTL(ind2)./normfac(ind2))+std(pCTL(ind2)./normfac(ind2))/sqrt(length(DSumm)-2)],'Linewidth',2)

        bar(2,mean(ACQ(ind2)./normfac(ind2)))
        plot([2 2],[mean(ACQ(ind2)./normfac(ind2))-std(ACQ(ind2)./normfac(ind2))/sqrt(length(DSumm)-2) mean(ACQ(ind2)./normfac(ind2))+std(ACQ(ind2)./normfac(ind2))/sqrt(length(DSumm)-2)],'Linewidth',2)


% pk20r49

% 1 % WN starts at dawn on 1.31
    % PRE % pk20r49/Exp1/acampoff_300108 % only the 1.31.08 data in this folder
    % DUR % /200muampon1.31.08
    % POST % /acampon13108-3
    
% 2 
    % PRE % pk20r49/Exp1/acampon13108-4 
    % DUR % /200muampon1.31.08
    % POST % acsf_ampon2
    
% 3 
    % PRE % pk20r49/Exp1/probein20408 
    % DUR % /200mu20508
    % POST % /ac020508
    
    num=3;
DSumm(num).fv{1}=fvPRE;
DSumm(num).fv{2}=fvDUR;
DSumm(num).fv{3}=fvPOST;
DSumm(num).pitch{3}=pitchPOST;
DSumm(num).pitch{2}=pitchDUR;
DSumm(num).pitch{1}=pitchPRE;   
    DSumm(num).bird='pk20r49';
    DSumm(num).window=190:210;
    figure;hold on;
    for i=1:length(DSumm(num).fv)
        for j=1:length(DSumm(num).fv{i})
            if DSumm(num).fv{i}(j).TRIG
                DSumm(num).istrig{i}(j)=1;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(window,j)),'.')
            else
                DSumm(num).istrig{i}(j)=0;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(window,j)),'r.')
            end
        end
    end
    
    % What is the CV?
        DSumm(num).CVpre=mean(jcstd(DSumm(num).pitch{1}(window,:)')./mean(DSumm(num).pitch{1}(window,:)'));
        DSumm(num).CVpost=mean(jcstd(DSumm(num).pitch{3}(window,:)')./mean(DSumm(num).pitch{3}(window,:)'));
    % What is the hit rate during the MUSC run?
        DSumm(num).HitRate=mean(DSumm(num).istrig{2});
    % What is the deltaFF?
        DSumm(num).deltaFF=mean(mean(DSumm(num).pitch{3}(window,:)')-mean(DSumm(num).pitch{1}(window,:)'));
        DSumm(num).mnFFpre=mean(mean(DSumm(num).pitch{1}(window,:)'));
        
        
% pu34

% 1 % 
    % PRE % ac716
    % DUR % /lid717
    % POST % /ac717
    
% 2 
    % PRE % 
    % DUR % 
    % POST % 
    
    num=5;
DSumm(num).coef=+1;
DSumm(num).fv{1}=fvPRE;
DSumm(num).fv{2}=fvDUR;
DSumm(num).fv{3}=fvPOST;
DSumm(num).pitch{3}=pitchPOST;
DSumm(num).pitch{2}=pitchDUR;
DSumm(num).pitch{1}=pitchPRE;   
    DSumm(num).bird='pu34';
    DSumm(num).window=604:644;
    figure;hold on;
    for i=1:length(DSumm(num).fv)
        for j=1:length(DSumm(num).fv{i})
            if DSumm(num).fv{i}(j).TRIG
                DSumm(num).istrig{i}(j)=1;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'.')
            else
                DSumm(num).istrig{i}(j)=0;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'r.')
            end
        end
    end
    
    % What is the CV?
        DSumm(num).CVpre=mean(jcstd(DSumm(num).pitch{1}(DSumm(num).window,:)')./mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));
        DSumm(num).CVpost=mean(jcstd(DSumm(num).pitch{3}(DSumm(num).window,:)')./mean(DSumm(num).pitch{3}(DSumm(num).window,:)'));
    % What is the hit rate during the MUSC run?
        DSumm(num).HitRate=mean(DSumm(num).istrig{2});
    % What is the deltaFF?
        DSumm(num).deltaFF=DSumm(num).coef*(mean(mean(DSumm(num).pitch{3}(DSumm(num).window,:)')-mean(DSumm(num).pitch{1}(DSumm(num).window,:)')));
        DSumm(num).mnFFpre=mean(mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));        
        
% bk20bk45

% 6 % -- A (7.15)
% 7 % -- B (7.15)
% 8 % -- A
% 9 % -- B
    % PRE % ac716
    % DUR % /lid717
    % POST % /ac717
    
% 2 
    % PRE % 
    % DUR % 
    % POST % 
    
    num=13;
DSumm(num).nt='b';
DSumm(num).coef=-1;
DSumm(num).fv{1}=fvPRE;
DSumm(num).fv{2}=fvDUR;
DSumm(num).fv{3}=fvPOST;
DSumm(num).pitch{3}=pitchPOST;
DSumm(num).pitch{2}=pitchDUR;
DSumm(num).pitch{1}=pitchPRE;   
    DSumm(num).bird='bk20bk45';
    DSumm(num).window=260:280;
    figure;hold on;
    for i=1:length(DSumm(num).fv)
        for j=1:length(DSumm(num).fv{i})
            if DSumm(num).fv{i}(j).TRIG
                DSumm(num).istrig{i}(j)=1;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'.')
            else
                DSumm(num).istrig{i}(j)=0;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'r.')
            end
        end
    end
    
    % What is the CV?
        DSumm(num).CVpre=mean(jcstd(DSumm(num).pitch{1}(DSumm(num).window,:)')./mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));
        DSumm(num).CVpost=mean(jcstd(DSumm(num).pitch{3}(DSumm(num).window,:)')./mean(DSumm(num).pitch{3}(DSumm(num).window,:)'));
    % What is the hit rate during the MUSC run?
        DSumm(num).HitRate=mean(DSumm(num).istrig{2});
    % What is the deltaFF?
        DSumm(num).deltaFF=DSumm(num).coef*(mean(mean(DSumm(num).pitch{3}(DSumm(num).window,:)')-mean(DSumm(num).pitch{1}(DSumm(num).window,:)')));
        DSumm(num).mnFFpre=mean(mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));        
        
% bk28w6

% 
    
    num=14;
DSumm(num).nt='b';
DSumm(num).coef=1;
DSumm(num).fv{1}=fvPRE;
DSumm(num).fv{2}=fvDUR;
DSumm(num).fv{3}=fvPOST;
DSumm(num).pitch{3}=pitchPOST;
DSumm(num).pitch{2}=pitchDUR;
DSumm(num).pitch{1}=pitchPRE;   
    DSumm(num).bird='bk28w6';
    DSumm(num).window=900:960;
    figure;hold on;
    for i=1:length(DSumm(num).fv)
        for j=1:length(DSumm(num).fv{i})
            if DSumm(num).fv{i}(j).TRIG
                DSumm(num).istrig{i}(j)=1;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'.')
            else
                DSumm(num).istrig{i}(j)=0;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'r.')
            end
        end
    end
    
    % What is the CV?
        DSumm(num).CVpre=mean(jcstd(DSumm(num).pitch{1}(DSumm(num).window,:)')./mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));
        DSumm(num).CVpost=mean(jcstd(DSumm(num).pitch{3}(DSumm(num).window,:)')./mean(DSumm(num).pitch{3}(DSumm(num).window,:)'));
    % What is the hit rate during the MUSC run?
        DSumm(num).HitRate=mean(DSumm(num).istrig{2});
    % What is the deltaFF?
        DSumm(num).deltaFF=DSumm(num).coef*(mean(mean(DSumm(num).pitch{3}(DSumm(num).window,:)')-mean(DSumm(num).pitch{1}(DSumm(num).window,:)')));
        DSumm(num).mnFFpre=mean(mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));        
       
% pk32bk28

% 
    
    num=15;
DSumm(num).nt='b';
DSumm(num).coef=1;
DSumm(num).fv{1}=fvPRE;
DSumm(num).fv{2}=fvDUR;
DSumm(num).fv{3}=fvPOST;
DSumm(num).pitch{3}=pitchPOST;
DSumm(num).pitch{2}=pitchDUR;
DSumm(num).pitch{1}=pitchPRE;   
    DSumm(num).bird='bk28w6';
    DSumm(num).window=900:960;
    figure;hold on;
    for i=1:length(DSumm(num).fv)
        for j=1:length(DSumm(num).fv{i})
            if DSumm(num).fv{i}(j).TRIG
                DSumm(num).istrig{i}(j)=1;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'.')
            else
                DSumm(num).istrig{i}(j)=0;
                plot(timing3(DSumm(num).fv{i}(j))-min(timing3(DSumm(num).fv{2})),...
                    mean(DSumm(num).pitch{i}(DSumm(num).window,j)),'r.')
            end
        end
    end
    
    % What is the CV?
        DSumm(num).CVpre=mean(jcstd(DSumm(num).pitch{1}(DSumm(num).window,:)')./mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));
        DSumm(num).CVpost=mean(jcstd(DSumm(num).pitch{3}(DSumm(num).window,:)')./mean(DSumm(num).pitch{3}(DSumm(num).window,:)'));
    % What is the hit rate during the MUSC run?
        DSumm(num).HitRate=mean(DSumm(num).istrig{2});
    % What is the deltaFF?
        DSumm(num).deltaFF=DSumm(num).coef*(mean(mean(DSumm(num).pitch{3}(DSumm(num).window,:)')-mean(DSumm(num).pitch{1}(DSumm(num).window,:)')));
        DSumm(num).mnFFpre=mean(mean(DSumm(num).pitch{1}(DSumm(num).window,:)'));        
               
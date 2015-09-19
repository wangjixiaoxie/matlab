% Earlier analysis % load /cardinal3/TrialbyTrial_Sept2010i.mat
% Most recent analysis
load /cardinal3/TrialbyTrial_Oct2010.mat

% Prediction window - 100 renditions where maximal learning occurs
% Learning measurement window - next 200 renditions (asymptote)

%%%% Plot individual data
%% Experiment #1 - bk80w28

        enum=1;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        figure;subplot(321);hold on;
        minpred=750;
        maxpred=850;
        learningwindow=[850:1050];
        ptwindow=[300:540];

        % Prediction and final learning blocks
        %plot(runningaverage(mean(Experiment(1).pitchWNon(ptwindow,[minpred:max(learningwindow)]))-mean(mean(Experiment(1).pitchBaseline(ptwindow,1:60))),20),'Linewidth',2)
        xlim([1 300])
        %plot([100 100],[1 300],'r')
        %
        indR=min(learningwindow)-1+find(std(ResidWN1(ptwindow,learningwindow))<0.02);
        indE1=indE(find(indE>minpred & indE<maxpred & std(ResidWN1(ptwindow,indE))<0.02));
        indH1=indH(find(indH>minpred & indH<maxpred & std(ResidWN1(ptwindow,indH))<0.02));
        %subplot(522);hold on;
        plot(mean(Experiment(1).pitchWNon(:,indR)')-mean(Experiment(1).pitchBaseline(:,1:60)'),'Linewidth',2)
        plot(mean(Experiment(1).pitchWNon(:,indE1)')-mean(Experiment(1).pitchBaseline(:,1:60)'),'k','Linewidth',2)
        plot(mean(Experiment(1).pitchWNon(:,indH1(2:end))')-mean(Experiment(1).pitchBaseline(:,1:60)'),'r','Linewidth',2)
        xlim([300 540])
        corrcoef(mean(Experiment(1).pitchWNon(ptwindow,indH1(2:end))')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)'),mean(Experiment(1).pitchWNon(ptwindow,indR)')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)'))
        corrcoef(mean(Experiment(1).pitchWNon(ptwindow,indE1)')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)'),mean(Experiment(1).pitchWNon(ptwindow,indR)')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)'))
        EscapePred(1).data=mean(Experiment(1).pitchWNon(ptwindow,indE1)')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)');
        HitPred(1).data=mean(Experiment(1).pitchWNon(ptwindow,indH1(2:end))')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)');
        Actual(1).data=mean(Experiment(1).pitchWNon(ptwindow,indR)')-mean(Experiment(1).pitchBaseline(ptwindow,1:60)');
        targ(1)=round(median(Experiment(1).targeting))-ptwindow(1);


%%%% Experiment # 2 - bk80w28
        ResidWN=ResidWN2;
        enum=2;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        subplot(322);hold on;
        minpred=600;
        maxpred=700;
        learningwindow=[700:900];
        ptwindow=[850:1100];

        % Prediction and final learning blocks
        %plot(runningaverage(mean(Experiment(2).pitchWNon(ptwindow,[minpred:max(learningwindow)]))-mean(mean(Experiment(2).pitchBaseline(ptwindow,:)))-100,20),'Linewidth',2)
        xlim([1 300])
        %plot([100 100],[1 300],'r')
        %
        indR=min(learningwindow)-1+find(std(ResidWN(ptwindow,learningwindow))<0.02);
        indE1=indE(find(indE>minpred & indE<maxpred & std(ResidWN(ptwindow,indE))<0.02));
        indH1=indH(find(indH>minpred & indH<maxpred & std(ResidWN(ptwindow,indH))<0.02));
        %subplot(524);hold on;
        plot(mean(Experiment(2).pitchWNon(ptwindow,indR)')-mean(Experiment(2).pitchBaseline(ptwindow,:)')-100,'Linewidth',2)
        plot(mean(Experiment(2).pitchWNon(ptwindow,indE1)')-mean(Experiment(2).pitchBaseline(ptwindow,:)')-100,'k','Linewidth',2)
        plot(mean(Experiment(2).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(2).pitchBaseline(ptwindow,:)')-100,'r','Linewidth',2)
        corrcoef(mean(Experiment(2).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(2).pitchBaseline(ptwindow,:)'),mean(Experiment(2).pitchWNon(ptwindow,indR)')-mean(Experiment(2).pitchBaseline(ptwindow,:)'))
        corrcoef(mean(Experiment(2).pitchWNon(ptwindow,indE1)')-mean(Experiment(2).pitchBaseline(ptwindow,:)'),mean(Experiment(2).pitchWNon(ptwindow,indR)')-mean(Experiment(2).pitchBaseline(ptwindow,:)'))

        EscapePred(2).data=mean(Experiment(2).pitchWNon(ptwindow,indE1)')-mean(Experiment(2).pitchBaseline(ptwindow,:)')-100;
        HitPred(2).data=mean(Experiment(2).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(2).pitchBaseline(ptwindow,:)')-100;
        Actual(2).data=mean(Experiment(2).pitchWNon(ptwindow,indR)')-mean(Experiment(2).pitchBaseline(ptwindow,:)')-100;
        targ(2)=round(median(Experiment(2).targeting))-ptwindow(1)+600;

%% Experiment #3 - bk80w28
        ResidWN=ResidWN3;
        enum=3;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        subplot(323);hold on;
        minpred=1;
        maxpred=100;
        learningwindow=[101:300];
        ptwindow=[300:515];

        % Prediction and final learning blocks
        %plot(runningaverage(mean(Experiment(3).pitchWNon(ptwindow,[minpred:max(learningwindow)]))-mean(mean(Experiment(3).pitchBaseline(ptwindow,[1:4 5:end]))),50),'Linewidth',2)
        %plot([100 100],[1 150],'r')
        %
        indR=min(learningwindow)-1+find(std(ResidWN(ptwindow,learningwindow))<0.02);
        indE1=indE(find(indE>minpred & indE<maxpred & std(ResidWN(ptwindow,indE))<0.02));
        indH1=indH(find(indH>minpred & indH<maxpred & std(ResidWN(ptwindow,indH))<0.02));
        % subplot(5,2,10);hold on;
        plot(mean(Experiment(3).pitchWNon(ptwindow,indR)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),'Linewidth',2)
        plot(mean(Experiment(3).pitchWNon(ptwindow,indE1)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),'k','Linewidth',2)
        plot(mean(Experiment(3).pitchWNon(ptwindow,indH1([1:5 7:end]))')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),'r','Linewidth',2)

        corrcoef(mean(Experiment(3).pitchWNon(ptwindow,indH1([1:5 7:end]))')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),mean(Experiment(3).pitchWNon(ptwindow,indR)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'))
        corrcoef(mean(Experiment(3).pitchWNon(ptwindow,indE1)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),mean(Experiment(3).pitchWNon(ptwindow,indR)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])'))

        EscapePred(3).data=mean(Experiment(3).pitchWNon(ptwindow,indE1)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])');
        HitPred(3).data=mean(Experiment(3).pitchWNon(ptwindow,indH1([1:5 7:end]))')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])');
        Actual(3).data=mean(Experiment(3).pitchWNon(ptwindow,indR)')-mean(Experiment(3).pitchBaseline(ptwindow,[1:4 6:35 37:end])');
        targ(3)=round(median(Experiment(3).targeting))-ptwindow(1);

%%%% Experiment # 4 - r87g80
        ResidWN=ResidWN4;
        enum=4;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        subplot(324);hold on;
        minpred=21;
        maxpred=120;
        learningwindow=[121:320];
        ptwindow=[250:450];

        % Prediction and final learning blocks
        %plot(runningaverage(mean(Experiment(4).pitchWNon(ptwindow,[minpred:max(learningwindow)]))-mean(mean(Experiment(4).pitchBaseline(ptwindow,:))),20),'Linewidth',2)
        %plot([100 100],[1 150],'r')
        %
        indR=min(learningwindow)-1+find(std(ResidWN(ptwindow,learningwindow))<0.2); % no outliers
        indE1=indE(find(indE>minpred & indE<maxpred & std(ResidWN(ptwindow,indE))<0.2));
        indH1=indH(find(indH>minpred & indH<maxpred & std(ResidWN(ptwindow,indH))<0.2));
        %subplot(526);hold on;
        plot(mean(Experiment(4).pitchWNon(ptwindow,indR)')-mean(Experiment(4).pitchBaseline(ptwindow,:)'),'Linewidth',2)
        plot(mean(Experiment(4).pitchWNon(ptwindow,indE1)')-mean(Experiment(4).pitchBaseline(ptwindow,:)'),'k','Linewidth',2)
        plot(mean(Experiment(4).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(4).pitchBaseline(ptwindow,:)'),'r','Linewidth',2)
        corrcoef(mean(Experiment(4).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(4).pitchBaseline(ptwindow,:)'),mean(Experiment(4).pitchWNon(ptwindow,indR)')-mean(Experiment(4).pitchBaseline(ptwindow,:)'))
        corrcoef(mean(Experiment(4).pitchWNon(ptwindow,indE1)')-mean(Experiment(4).pitchBaseline(ptwindow,:)'),mean(Experiment(4).pitchWNon(ptwindow,indR)')-mean(Experiment(4).pitchBaseline(ptwindow,:)'))
 
        EscapePred(4).data=mean(Experiment(4).pitchWNon(ptwindow,indE1)')-mean(Experiment(4).pitchBaseline(ptwindow,:)');
        HitPred(4).data=mean(Experiment(4).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(4).pitchBaseline(ptwindow,:)');
        Actual(4).data=mean(Experiment(4).pitchWNon(ptwindow,indR)')-mean(Experiment(4).pitchBaseline(ptwindow,:)');
        targ(4)=0;%round(median(Experiment(4).targeting))-ptwindow(1);

%%%% Experiment # 5 - bk80w28
        ResidWN=ResidWN5;
        enum=5;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        subplot(325);hold on;
        minpred=241;
        maxpred=341;
        learningwindow=[341:541];
        ptwindow=[300:550];

        % Prediction and final learning blocks
        %plot(runningaverage(mean(Experiment(5).pitchWNon(ptwindow,[minpred:max(learningwindow)]))-mean(mean(Experiment(5).pitchBaseline(ptwindow,[1:4 5:end]))),10),'Linewidth',2)
        %plot([100 100],[1 150],'r')
        %
        indR=min(learningwindow)-1+find(std(ResidWN(ptwindow,learningwindow))<0.2); % no outliers
        indE1=indE(find(indE>minpred & indE<maxpred & std(ResidWN(ptwindow,indE))<0.2));
        indH1=indH(find(indH>minpred & indH<maxpred & std(ResidWN(ptwindow,indH))<0.2));
        %subplot(528);hold on;
        plot(median(Experiment(5).pitchWNon(ptwindow,indR)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),'Linewidth',2)
        plot(mean(Experiment(5).pitchWNon(ptwindow,indE1)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),'k','Linewidth',2)
        plot(mean(Experiment(5).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),'r','Linewidth',2)
        corrcoef(mean(Experiment(5).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),median(Experiment(5).pitchWNon(ptwindow,indR)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'))
        corrcoef(mean(Experiment(5).pitchWNon(ptwindow,indE1)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'),median(Experiment(5).pitchWNon(ptwindow,indR)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])'))

        EscapePred(5).data=mean(Experiment(5).pitchWNon(ptwindow,indE1)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])');
        HitPred(5).data=mean(Experiment(5).pitchWNon(ptwindow,indH1(1:end))')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])');
        Actual(5).data=mean(Experiment(5).pitchWNon(ptwindow,indR)')-mean(Experiment(5).pitchBaseline(ptwindow,[1:4 6:35 37:end])');
        targ(5)=round(median(Experiment(5).targeting))-ptwindow(1);

% Experiment #6 - bk91w60
        ResidWN=ResidWN6;
        enum=6;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        subplot(326);hold on;
        minpred=401;
        maxpred=501;
        learningwindow=[500:700];
        ptwindow=[255:345];
        % Prediction and final learning blocks
        %plot(runningaverage(mean(Experiment(6).pitchWNon(ptwindow,[minpred:max(learningwindow)]))-mean(mean(Experiment(6).pitchBaseline(ptwindow,:))),20),'Linewidth',2)
        %plot([100 100],[1 150],'r')
        %
        indR=min(learningwindow)-1+find(std(ResidWN(ptwindow,learningwindow))<0.02); % no outliers
        indE1=indE(find(indE>minpred & indE<maxpred & std(ResidWN(ptwindow,indE))<0.02));
        indH1=indH(find(indH>minpred & indH<maxpred & std(ResidWN(ptwindow,indH))<0.02));
        plot(mean(Experiment(6).pitchWNon(ptwindow,indR)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'),'Linewidth',2)
        plot(mean(Experiment(6).pitchWNon(ptwindow,indE1)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'),'k','Linewidth',2)
        plot(mean(Experiment(6).pitchWNon(ptwindow,indH1)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'),'r','Linewidth',2)
        corrcoef(mean(Experiment(6).pitchWNon(ptwindow,indH1)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'),mean(Experiment(6).pitchWNon(ptwindow,indR)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'))
        corrcoef(mean(Experiment(6).pitchWNon(ptwindow,indE1)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'),mean(Experiment(6).pitchWNon(ptwindow,indR)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])'))
      
        EscapePred(6).data=mean(Experiment(6).pitchWNon(ptwindow,indE1)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])');
        HitPred(6).data=mean(Experiment(6).pitchWNon(ptwindow,indH1)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])');
        Actual(6).data=mean(Experiment(6).pitchWNon(ptwindow,indR)')-mean(Experiment(6).pitchBaseline(ptwindow,[1:38 41:end])');
        targ(6)=round(median(Experiment(6).targeting))-ptwindow(1);


                
                
                

            
            
            

%%%%%%%% Summary quantification
            % Normalize to maximum
            figure;hold on;
            for i=1:6
                subplot(2,3,i);hold on;
                plot(Actual(i).data/max(Actual(i).data))
                plot(HitPred(i).data/max(HitPred(i).data),'r')
                plot(EscapePred(i).data/max(EscapePred(i).data),'k')
                distH(i).data=(abs(Actual(i).data/max(Actual(i).data)-HitPred(i).data/max(HitPred(i).data)));
                 distE(i).data=(abs(Actual(i).data/max(Actual(i).data)-EscapePred(i).data/max(EscapePred(i).data)));
                 dHm(i)=mean(distH(i).data);
                 dEm(i)=mean(distE(i).data);
                 ylim([0 1])
            end

            figure;hold on;
            subplot(2,3,1);hold on;
            i=3;
                plot(Actual(i).data/max(Actual(i).data),'Linewidth',2)
                plot(HitPred(i).data/max(HitPred(i).data),'r','Linewidth',2)
                plot(EscapePred(i).data/max(EscapePred(i).data),'k','Linewidth',2)
                ylim([0 1])
            subplot(2,3,2);hold on;
            i=5;
                plot(Actual(i).data/max(Actual(i).data),'Linewidth',2)
                plot(HitPred(i).data/max(HitPred(i).data),'r','Linewidth',2)
                plot(EscapePred(i).data/max(EscapePred(i).data),'k','Linewidth',2)
                ylim([0 1])

            subplot(2,3,3);hold on;
            plot(ones(1,6),dHm,'.','Markersize',15,'Color','r')
            plot(2*ones(1,6),dEm,'.','Markersize',15,'Color','k')
            plot([ones(1,6);2*ones(1,6)],[dHm;dEm],'k')
            xlim([0 3])


            [h,p]=ttest(dHm,dEm)

            for i=1:5
                H=corrcoef(Actual(i).data,HitPred(i).data);
                E=corrcoef(Actual(i).data,EscapePred(i).data);
                 cHm(i)=H(2);
                 cEm(i)=E(2);
            end
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

































figure;hold on;
for i=1:5
    plot(distH(i).data-distE(i).data,'r')
end
% Align to targeting time and get standard error!
alignA=zeros(5,500);
for i=1:5
    tn=targ(i);
    alignA(i,501-tn:500)=distH(i).data(1:tn)-distE(i).data(1:tn);
    nextchunk=distH(i).data(tn+1:end)-distE(i).data(tn+1:end);
    alignA(i,501:500+length(nextchunk))=nextchunk;
end
win=20;
  figure;plot(runningaverage(mean(alignA(:,440:670)),win))
hold on;plot(runningaverage(mean(alignA(:,440:670)),win)-runningaverage(std(alignA(:,440:670)),win)/sqrt(5),'r')
hold on;plot([0 200],[0 0],'k')      
       
alignH=zeros(5,700);
for i=1:5
    tn=targ(i);
    alignH(i,501-tn:500)=HitPred(i).data(1:tn)/max(HitPred(i).data(1:tn));
    nextchunk=HitPred(i).data(tn+1:end)/max(HitPred(i).data(tn+1:end));
    alignH(i,501:500+length(nextchunk))=nextchunk;
end
alignE=zeros(5,700);
for i=1:5
    tn=targ(i);
    alignE(i,501-tn:500)=EscapePred(i).data(1:tn)/max(EscapePred(i).data(1:tn));
    nextchunk=EscapePred(i).data(tn+1:end)/max(EscapePred(i).data(tn+1:end));
    alignE(i,501:500+length(nextchunk))=nextchunk;
end
alignAC=zeros(5,700);
for i=1:5
    tn=targ(i);
    alignAC(i,501-tn:500)=Actual(i).data(1:tn)/max(Actual(i).data(1:tn));
    nextchunk=Actual(i).data(tn+1:end)/max(Actual(i).data(tn+1:end));
    alignAC(i,501:500+length(nextchunk))=nextchunk;
end
for i=1:700
    MalignH(i)=mean(alignH(find(alignH(:,i)~=0),i));
    MalignE(i)=mean(alignE(find(alignE(:,i)~=0),i));
    MalignAC(i)=mean(alignAC(find(alignAC(:,i)~=0),i));
end




figure;hold on;
for i=1:5
subplot(2,3,i);hold on;
plot(Actual(i).data/mean(Actual(i).data))
plot(HitPred(i).data/mean(HitPred(i).data),'r')
plot(EscapePred(i).data/mean(EscapePred(i).data),'k')
distH(i).data=(abs(Actual(i).data/mean(Actual(i).data)-HitPred(i).data/mean(HitPred(i).data)));
distE(i).data=(abs(Actual(i).data/mean(Actual(i).data)-EscapePred(i).data/mean(EscapePred(i).data)));
dHm(i)=mean(distH(i).data);
dEm(i)=mean(distE(i).data);
end
win=20;
figure;plot(runningaverage(mean(alignB),win))
hold on;plot(runningaverage(mean(alignB),win)-runningaverage(std(alignB),win)/sqrt(5),'r')
hold on;plot([0 200],[0 0],'k')

        
        
        
        
        
        
figure;plot(mean(ResidWN6(300:450,indE(1:15))'),'r')
hold on;plot(mean(ResidWN6(300:450,300:end)'),'k')
hold on;plot(mean(ResidWN6(300:450,indH(1:15))'),'b')
% significant for experiment #1 - only during learning portion ~(750:950)
        figure;plot(runningaverage(mean(Experiment(1).pitchWNon(300:550,650:end)),20))
        enum=1;
                indE=find(Experiment(enum).indEscapeAbove==1);
                indH=find(Experiment(enum).indHitAbove==1);
        clear cH cE
        for i=1:1163
            indH1=indH(find(indH>i & indH<i+100 & std(ResidWN1(300:550,indH))<0.02));
            indR=i+100+find(std(ResidWN1(300:550,i+100:i+200))<0.02);
            c1=corrcoef(mean(ResidWN1(300:550,indH1)'),mean(ResidWN1(300:550,indR)'));
            indE1=indE(find(indE>i & indE<i+100 & std(ResidWN1(300:600,indE))<0.02));
            indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
            c2=corrcoef(mean(ResidWN1(300:550,indE2)'),mean(ResidWN1(300:550,indR)'));
            cH(i)=c1(2);
            cE(i)=c2(2);
        end
        figure;plot(cE(650:end));hold on;plot(cH(650:end),'r')


% significant for experiment #2 - throughout
        figure;plot(runningaverage(mean(Experiment(2).pitchWNon(900:1100,1:end)),20))

                ResidWN=ResidWN2;
                enum=2;
                indE=find(Experiment(enum).indEscapeAbove==1);
                indH=find(Experiment(enum).indHitAbove==1);
        clear cH cE
        for i=1:934
            indH1=indH(find(indH>i & indH<i+100 & std(ResidWN(900:1100,indH))<0.025));
            indR=i+100+find(std(ResidWN(900:1100,i+101:i+200))<0.025);
            c1=corrcoef(mean(ResidWN(900:1100,indH1)'),mean(ResidWN(900:1100,indR)'));
            indE1=indE(find(indE>i & indE<i+100 & std(ResidWN(900:1100,indE))<0.025));
            indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
            c2=corrcoef(mean(ResidWN(900:1100,indE2)'),mean(ResidWN(900:1100,indR)'));
            cH(i)=c1(2);
            cE(i)=c2(2);
        end
        figure;plot(cE);hold on;plot(cH,'r')

% significant for experiment #6 - only in learning portion
        figure;plot(runningaverage(mean(Experiment(6).pitchWNon(300:450,1:end)),20))
                ResidWN=ResidWN6;
                enum=6;
                indE=find(Experiment(enum).indEscapeAbove==1);
                indH=find(Experiment(enum).indHitAbove==1);
        clear cH cE
        for i=1:184
            indH1=indH(find(indH>i & indH<i+100 & std(ResidWN6(300:450,indH))<0.5));
            indR=i+100+find(std(ResidWN6(300:450,i+101:i+200))<0.5);
            c1=corrcoef(mean(ResidWN6(300:450,indH1)'),mean(ResidWN6(300:450,indR)'));
            indE1=indE(find(indE>i & indE<i+100 & std(ResidWN6(300:450,indE))<0.5));
            indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
            c2=corrcoef(mean(ResidWN6(300:450,indE2)'),mean(ResidWN6(300:450,indR)'));
            cH(i)=c1(2);
            cE(i)=c2(2);
        end
        figure;plot(cE);hold on;plot(cH,'r')

%%%%%%%%%%%%
        figure;plot(runningaverage(mean(Experiment(5).pitchWNon(300:550,1:end)),100))
% ~600-900 is the only epoch of learning where the notch is decent
                ResidWN=ResidWN5;
                enum=5;
                indE=find(Experiment(enum).indEscapeAbove==1);
                indH=find(Experiment(enum).indHitAbove==1);
        clear cH cE
        for i=1:930 % up until notch starts being shitty again
            indH1=indH(find(indH>i & indH<i+100 & std(ResidWN(300:550,indH))<0.025));
            indR=i+100+find(std(ResidWN(300:550,i+101:i+200))<0.025);
            c1=corrcoef(mean(ResidWN(300:550,indH1)'),mean(ResidWN(300:550,indR)'));
            indE1=indE(find(indE>i & indE<i+100 & std(ResidWN(300:550,indE))<0.025));
            indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
            c2=corrcoef(mean(ResidWN(300:550,indE2)'),mean(ResidWN(300:550,indR)'));
            cH(i)=c1(2);
            cE(i)=c2(2);
        end
        figure;plot(cE);hold on;plot(cH,'r')




        figure;plot(runningaverage(mean(Experiment(7).pitchWNon(300:450,1:end)),100))
% ~600-900 is the only epoch of learning where the notch is decent
                ResidWN=ResidWN7;
                enum=7;
                indE=find(Experiment(enum).indEscapeAbove==1);
                indH=find(Experiment(enum).indHitAbove==1);
        clear cH cE
        for i=1:930 % up until notch starts being shitty again
            indH1=indH(find(indH>i & indH<i+100 & std(ResidWN(250:450,indH))<0.025));
            indR=i+100+find(std(ResidWN(250:450,i+101:i+200))<0.025);
            c1=corrcoef(mean(ResidWN(250:450,indH1)'),mean(ResidWN(250:450,indR)'));
            indE1=indE(find(indE>i & indE<i+100 & std(ResidWN(250:450,indE))<0.025));
            indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
            c2=corrcoef(mean(ResidWN(250:450,indE2)'),mean(ResidWN(250:450,indR)'));
            cH(i)=c1(2);
            cE(i)=c2(2);
        end
        figure;plot(cE);hold on;plot(cH,'r')





% 
figure;hold on;
subplot(331);
plot(runningaverage(mean(Experiment(3).pitchWNon(310:550,250:end)),100),'Linewidth',2)
xlim([0 210])
subplot(334);hold on;
plot(runningaverage(cccE3,10),'k','Linewidth',2)
plot(runningaverage(cccH3,10),'r','Linewidth',2)
subplot(332);
plot(runningaverage(mean(Experiment(6).pitchWNon(300:450,1:200)),100),'Linewidth',2)
subplot(335);hold on;
plot(runningaverage(cccE6,10),'k','Linewidth',2)
plot(runningaverage(cccH6,10),'r','Linewidth',2)
ylim([0 0.2])
subplot(333);
plot(runningaverage(mean(Experiment(7).pitchWNon(300:450,130:end)),100),'Linewidth',2)
xlim([0 130])
subplot(336);hold on;
plot(runningaverage(cccE7,10),'k','Linewidth',2)
plot(runningaverage(cccH7,10),'r','Linewidth',2)
ylim([0 0.2])
subplot(337);hold on;
        plot(runningaverage([m3;m6;m7],10))
        plot(runningaverage([m3;m6;m7],10)-runningaverage(std([m3;m6;m7]),10)/sqrt(3),'Linewidth',2)
        plot(runningaverage([m3;m6;m7],10)+runningaverage(std([m3;m6;m7]),10)/sqrt(3),'Linewidth',2)        
        plot([0 50],[0 0],'k')
        
        
        
% All examples are highly significant - group data are consistent
        figure;hold on;
        plot(cccE3-cccH3,'k','Linewidth',2)
        plot(cccE6-cccH6,'b','Linewidth',2)
        plot(cccE7-cccH7,'r','Linewidth',2)
        m3=cccE3-cccH3;
        m6=cccE6-cccH6;
        m7=cccE7-cccH7;
        figure;hold on;
        plot(runningaverage([m3;m6;m7],10))
        plot(runningaverage([m3;m6;m7],10)-runningaverage(std([m3;m6;m7]),10)/sqrt(3))
        plot(runningaverage([m3;m6;m7],10)+runningaverage(std([m3;m6;m7]),10)/sqrt(3))        
        plot([0 50],[0 0],'k')
        figure;hold on;
        plot(mean([m3;m6;m7]))
        plot(mean([m3;m6;m7])-std([m3;m6;m7])/sqrt(3))
        plot(mean([m3;m6;m7])+std([m3;m6;m7])/sqrt(3))
        plot([0 50],[0 0],'k')
% Figure for experiment 3
        enum=3;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
clear t
for i=1:202
    t(i,:)=1/8:1/8:241/8;
end
figure;hold on;
subplot(321);hold on; % When does learning occur
plot(runningaverage(mean(Experiment(3).pitchWNon(310:550,:)),100),'Linewidth',2)
subplot(322);hold on;
plot(t',Experiment(3).pitchWNon(310:550,indH),'r')
plot(t(1:191,:)',Experiment(3).pitchWNon(310:550,indE),'k')
    xlim([1/8 241/8])
ylim([6600 7500])   
    



% Experiment 1
figure;hold on;
% s.d. < 0.015 gives you non-outliers

indH1=indH(find(indH>800 & indH<900 & std(ResidWN1(300:500,indH))<0.015));
corrcoef(median(ResidWN1(300:500,indH1)'),median(ResidWN1(300:500,900:1000)'))
indE1=indE(find(indE>800 & indE<900 & std(ResidWN1(300:500,indE))<0.015));
indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
corrcoef(median(ResidWN1(300:500,indE2)'),median(ResidWN1(300:500,900:1000)'))

indH1=indH(find(indH>300 & indH<400 & std(ResidWN1(300:500,indH))<0.015));
corrcoef(median(ResidWN1(300:500,indH1)'),median(ResidWN1(300:500,400:500)'))
indE1=indE(find(indE>300 & indE<400 & std(ResidWN1(300:500,indE))<0.015));
indE2=indE1;%(ceil(rand(1,12)*length(indH1)));
corrcoef(median(ResidWN1(300:500,indE2)'),median(ResidWN1(300:500,400:500)'))







% In the chunk with learning
% Learning epoch
        ResidWN=ResidWN6(:,1:200);
        enum=1;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        indH=indH(find(indH>0));
        indE=indE(find(indE>0));
        clear cccH cccE ccH ccE
        window1=[310:450];
        for kk=1:50 % For each time in the future
            for i=1:length(indH(indH+50<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i,kk)=a(2);
            end
            for i=1:length(indE(indE+50<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i,kk)=b(2);
            end
            cccH(kk)=mean(ccH(:,kk));
            cccE(kk)=mean(ccE(:,kk));
        end
        mean(cccE-cccH) % 0.784 --- p<1e-20
%%%%%%%%%
    subplot(323);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:241/8,mean(Experiment(3).pitchWNon(310:550,indE(find(indE>250 & indE<500)))'),'k','Linewidth',2)
    plot(1/8:1/8:241/8,mean(Experiment(3).pitchWNon(310:550,indH(find(indH>250 & indH<500)))'),'r','Linewidth',2)
    plot(1/8:1/8:241/8,mean(Experiment(3).pitchWNon(310:550,350:550)'),'b','Linewidth',2)
    xlim([1/8 241/8])
    ylim([7000 7100])
    subplot(325);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,3),'k','Linewidth',2);
    plot(runningaverage(cccH,3),'r','Linewidth',2);
    ylim([-0.05 0.2])
%     subplot(336);hold on; % Which H's and E's are responsible for this difference?
%     plot(runningaverage(mean(ccE'),10),'k','Linewidth',2);
%     plot(runningaverage(mean(ccH'),10),'r','Linewidth',2);
% In the chunk without learning
% Non-learning epoch
        ResidWN=ResidWN3(:,1:246);
        enum=3;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[310:550];
        for kk=1:50 % For each time in the future
            for i=1:length(indH(indH+50<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i,kk)=a(2);
            end
            for i=1:length(indE(indE+50<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i,kk)=b(2);
            end
            cccH(kk)=mean(ccH(:,kk));
            cccE(kk)=mean(ccE(:,kk));
        end
        mean(cccE-cccH) % -0.0265 --- p=0.0048
%%%%%%%%%%%%%
    subplot(324);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:251/8,mean(Experiment(3).pitchWNon(300:550,indE(find(indE>1 & indE<150)))'),'k','Linewidth',2)
    plot(1/8:1/8:251/8,mean(Experiment(3).pitchWNon(300:550,indH(find(indH>1 & indH<150)))'),'r','Linewidth',2)
    plot(1/8:1/8:251/8,mean(Experiment(3).pitchWNon(300:550,50:200)'),'b','Linewidth',2)
    xlim([1/8 241/8])
    ylim([6900 7100])
    subplot(326);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,10),'k','Linewidth',2);
    plot(runningaverage(cccH,10),'r','Linewidth',2);
    ylim([-0.05 0.2])
%     subplot(339);hold on; % Which H's and E's are responsible for this difference?
%     plot(runningaverage(mean(ccE'),10),'k','Linewidth',2);
%     plot(runningaverage(mean(ccH'),10),'r','Linewidth',2);


% Figure for Experiment 6
        enum=6;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
clear t
for i=1:130
    t(i,:)=1/8:1/8:151/8;
end
figure;hold on;
subplot(321);hold on; % When does learning occur
plot(runningaverage(mean(Experiment(6).pitchWNon(310:450,:)),100),'Linewidth',2)
xlim([0 280])
subplot(322);hold on;
plot(t(1:128,:)',Experiment(6).pitchWNon(300:450,indH),'r')
plot(t(1:130,:)',Experiment(6).pitchWNon(300:450,indE),'k')
    xlim([1/8 151/8])
ylim([7500 8500])   
%%%%%%%%%
        ResidWN=ResidWN4(:,1:500);
        enum=4;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:50 % For each time in the future
            for i=1:length(indH(indH+50<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+50<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % 0.0585 --- p<1e-6
%%%%%%%%
    subplot(323);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,indE(find(indE>1 & indE<150)))'),'k','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,indH(find(indH>1 & indH<150)))'),'r','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,50:200)'),'b','Linewidth',2)
    xlim([1/8 151/8])
    %ylim([7000 7100])
    subplot(325);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,10),'k','Linewidth',2);
    plot(runningaverage(cccH,10),'r','Linewidth',2);
    ylim([0 0.2])
%%%%%%%%%%
        ResidWN=ResidWN6(:,201:end);
        expnum=6;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:50 % For each time in the future
            for i=1:length(indH(indH+50<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+50<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % 0.0324 --- p<0.01
%%%%%%%%%%%%%
    subplot(324);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,indE(find(indE>200 & indE<300)))'),'k','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,indH(find(indH>200 & indH<300)))'),'r','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,250:350)'),'b','Linewidth',2)
    xlim([1/8 151/8])
    %ylim([7000 7100])
    subplot(326);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,10),'k','Linewidth',2);
    plot(runningaverage(cccH,10),'r','Linewidth',2);
    ylim([0 0.2])


% Figure for Experiment 7
        enum=7;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
clear t
for i=1:130
    t(i,:)=1/8:1/8:151/8;
end
figure;hold on;
subplot(321);hold on; % When does learning occur
plot(runningaverage(mean(Experiment(7).pitchWNon(310:450,:)),100),'Linewidth',2)
xlim([0 280])
subplot(322);hold on;
plot(t(1:107,:)',Experiment(7).pitchWNon(300:450,indH),'r')
plot(t(1:100,:)',Experiment(7).pitchWNon(300:450,indE),'k')
    xlim([1/8 151/8])
ylim([7600 8600])   


% Learning epoch
        ResidWN=ResidWN7(:,130:end);
        enum=7;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:50 % For each time in the future
            for i=1:length(indH(indH+50<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+50<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % 0.0285 --- p<0.002
 %%%%%%%%%%%%
     subplot(323);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:151/8,mean(Experiment(7).pitchWNon(300:450,indE(find(indE>181 & indE<300)))'),'k','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(7).pitchWNon(300:450,indH(find(indH>181 & indH<300)))'),'r','Linewidth',2) % 181 because 180 is outlier
    plot(1/8:1/8:151/8,mean(Experiment(7).pitchWNon(300:450,200:350)'),'b','Linewidth',2)
    xlim([1/8 151/8])
    %ylim([7000 7100])
    subplot(325);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,10),'k','Linewidth',2);
    plot(runningaverage(cccH,10),'r','Linewidth',2);
    ylim([0 0.2])

% Non-learning epoch
        ResidWN=ResidWN7(:,1:130);
        enum=7;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:50 % For each time in the future
            for i=1:length(indH(indH+50<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+50<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % -0.033 --- p=0.05
%%%%%%%%%%%
    subplot(324);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,indE(find(indE>1 & indE<130)))'),'k','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,indH(find(indH>1 & indH<130)))'),'r','Linewidth',2)
    plot(1/8:1/8:151/8,mean(Experiment(6).pitchWNon(300:450,50:180)'),'b','Linewidth',2)
    xlim([1/8 151/8])
    %ylim([7000 7100])
    subplot(326);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,10),'k','Linewidth',2);
    plot(runningaverage(cccH,10),'r','Linewidth',2);
    ylim([0 0.2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% When learning is not happening, escapes and hits equally predict the future...
%%%% When learning happens, escapes above predict the future...
%%%% Even though mean hits and escapes are identical

%%%%%% To what extent does the structure of each escape/hit predict the nth future
%%%%%% performance?

% Experiment #3
% Learning epoch
        ResidWN=ResidWN3(:,250:end);
        enum=3;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[310:550];
        for kk=1:100 % For each time in the future
            for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i,kk)=a(2);
            end
            for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i,kk)=b(2);
            end
            cccH(kk)=mean(ccH(:,kk));
            cccE(kk)=mean(ccE(:,kk));
        end
        mean(cccE-cccH) % 0.784 --- p<1e-20
% Non-learning epoch
        ResidWN=ResidWN3(:,1:246);
        enum=3;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[310:550];
        for kk=1:100 % For each time in the future
            for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i,kk)=a(2);
            end
            for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i,kk)=b(2);
            end
            cccH(kk)=mean(ccH(:,kk));
            cccE(kk)=mean(ccE(:,kk));
        end
        mean(cccE-cccH) % -0.0265 --- p=0.0048
% Experiment #6 
% Learning epoch
        ResidWN=ResidWN6(:,1:200);
        enum=6;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:100 % For each time in the future
            for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % 0.0585 --- p<1e-6
% Non-learning epoch
        ResidWN=ResidWN6(:,151:300);
        expnum=6;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:100 % For each time in the future
            for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % -0.0349 --- p<0.03

% Experiment #7
% Learning epoch
        ResidWN=ResidWN7(:,130:end);
        enum=7;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:100 % For each time in the future
            for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % 0.0285 --- p<0.002
        
% Non-learning epoch
        ResidWN=ResidWN7(:,1:130);
        enum=7;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1);
        clear cccH cccE ccH ccE
        window1=[300:450];
        for kk=1:100 % For each time in the future
            for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
                a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                ccH(i)=a(2);
            end
            for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
                 b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                ccE(i)=b(2);
            end
            cccH(kk)=mean(ccH);
            cccE(kk)=mean(ccE);
        end
        mean(cccE-cccH) % -0.033 --- p=0.05




%%%%%%%%%%%%
%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%



indH=find(Experiment(enum).indHitAbove==1);
indE=find(Experiment(enum).indEscapeAbove==1);
indBoth=[indH indE];
window1=[300:450];
clear cE cH CDiff
window2=100:-1:1;
for i=window2(1)+1:size(ResidWN,2)
    mE=mean(ResidWN(window1,find(Experiment(enum).indEscapeAbove(i-window2)==1))');
    mH=mean(ResidWN(window1,find(Experiment(enum).indHitAbove(i-window2)==1))');
    a=corrcoef(ResidWN(window1,i),mE);
    cE(i-window2(1))=a(2);
    b=corrcoef(ResidWN(window1,i),mH);
    cH(i-window2(1))=b(2);    
end
CDiff=mean(cE-cH);% Predict greater than zero

figure;plot(runningaverage(cE,100))
hold on;plot(runningaverage(cH,100),'r')


% Calculate residuals
enum=3;

clear ResidWN
for i=1:size(Experiment(enum).pitchWNon,2)
    ResidWN(:,i)=(Experiment(enum).pitchWNon(:,i)./mean(Experiment(enum).pitchBaseline')')-1;
end
enum=3;
enum=1;
indE=find(Experiment(enum).indEscapeAbove==1);
indH=find(Experiment(enum).indHitAbove==1);





%%%%%%% To what extent does the structure of the nth escape/hit predict
%%%%%%% future performance?
clear cccH cccE ccH ccE
window1=[300:450];

    for i=1:length(indH(indH+100<size(ResidWN,2))) % For each hit
        for kk=1:100 % For each hit
            a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
            ccH(kk)=a(2);
        end
        cccH(i)=mean(ccH);
    end

    for i=1:length(indE(indE+100<size(ResidWN,2))) % For each escape
        for kk=1:100 % For each hit
         b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
        ccE(kk)=b(2);   
        end
        cccE(i)=mean(ccE);
    end



[h,p]=ttest(cccE-cccH) % for #6 it is 1e-10



figure;plot(runningaverage(cccE,20))
hold on;plot(runningaverage(cccH,20),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% When learning is not happening, escapes and hits equally predict the future...
%%%% When learning happens, escapes above predict the future...
%%%% Even though mean hits and escapes are identical

figure;hold on;
subplot(321);hold on;
plot(Experiment(6).pitchWNon(:,Experiment(6).indEscapeAbove),'b')
plot(Experiment(6).pitchWNon(:,Experiment(6).indHitAbove),'r')
xlim([200 550])
subplot(322);hold on;
plot(Experiment(7).pitchWNon(:,Experiment(7).indEscapeAbove),'b')
plot(Experiment(7).pitchWNon(:,Experiment(7).indHitAbove),'r')
xlim([200 550])




        enum=1;
        ResidWN=ResidWN3;
        indE=find(Experiment(enum).indEscapeAbove==1);
        indH=find(Experiment(enum).indHitAbove==1); 
        indE=indE(find(indE>740 & indE<980));
        indH=indH(find(indH>740 & indH<980));        
        clear cccH cccE ccH ccE
        window1=[300:500];
        for kk=1:50 % For each time in the future
            countH=0;
            for i=1:length(indH(indH+50<980)) % For each hit
                if indsOUTLIER~=indH(i)+kk
                    countH=countH+1;
                    a=corrcoef(ResidWN(window1,indH(i)),ResidWN(window1,indH(i)+kk));
                    ccH(countH,kk)=a(2);
                end
            end
            countE=0;
            for i=1:length(indE(indE+50<980)) % For each escape
                if indsOUTLIER~=indE(i)+kk
                    countE=countE+1;
                    b=corrcoef(ResidWN(window1,indE(i)),ResidWN(window1,indE(i)+kk));
                    ccE(countE,kk)=b(2);
                end
            end
            cccH(kk)=mean(ccH(:,kk));
            cccE(kk)=mean(ccE(:,kk));
        end
        mean(cccE-cccH) % 0.784 --- p<1e-20












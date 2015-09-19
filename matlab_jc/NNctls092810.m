load /cardinal3/TrialbyTrial_Sept2010e.mat


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
    
% In the chunk with learning
% Learning epoch
        ResidWN=ResidWN3(:,250:end);
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
        mean(cccE-cccH) % 0.784 --- p<1e-20
%%%%%%%%%
    subplot(323);hold on; % What does mean H vs E look like there?
    plot(1/8:1/8:241/8,mean(Experiment(3).pitchWNon(310:550,indE(find(indE>250 & indE<500)))'),'k','Linewidth',2)
    plot(1/8:1/8:241/8,mean(Experiment(3).pitchWNon(310:550,indH(find(indH>250 & indH<500)))'),'r','Linewidth',2)
    plot(1/8:1/8:241/8,mean(Experiment(3).pitchWNon(310:550,350:550)'),'b','Linewidth',2)
    xlim([1/8 241/8])
    ylim([7000 7100])
    subplot(325);hold on; % How is future learning predicted differently by H vs E?
    plot(runningaverage(cccE,10),'k','Linewidth',2);
    plot(runningaverage(cccH,10),'r','Linewidth',2);
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
        ResidWN=ResidWN6(:,1:200);
        enum=6;
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


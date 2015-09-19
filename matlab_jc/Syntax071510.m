
% Syntax analysis 7-15-10


% r39g39
load /cardinal3/r39g9/SummaryPM.mat

for i=1:length(r39g39syn)
    r39g39syn(i).LClearn=r39g39syn(i).labels(find(r39g39syn(i).labels=='a' | r39g39syn(i).labels=='b' | r39g39syn(i).labels=='c' | r39g39syn(i).labels=='i'));
    r39g39syn(i).LCtargets=(r39g39syn(i).LClearn=='a' | r39g39syn(i).LClearn=='c');
end
% baseline - acsf
pACbase=mean([r39g39syn(1).LCtargets;r39g39syn(3).LCtargets]) % 0.7638
pAPVbase=mean(r39g39syn(2).LCtargets) % 0.8113
tAPVbase=mean(timing3([r39g39syn(1).fvalsA r39g39syn(1).fvalsB r39g39syn(1).fvalsC r39g39syn(1).fvalsI])); % 0.7674
pAPVlearn1=mean(r39g39syn(5).LCtargets) % 0.6974
pAPVlearn2=mean(r39g39syn(7).LCtargets) % 0.5517
pAPVlearn3=mean(r39g39syn(9).LCtargets) % 0.4571
tAPVlearn1=mean(timing3([r39g39syn(5).fvalsA r39g39syn(5).fvalsB r39g39syn(5).fvalsC r39g39syn(5).fvalsI])); % 0.7674
tAPVlearn2=mean(timing3([r39g39syn(7).fvalsA r39g39syn(7).fvalsB r39g39syn(7).fvalsC r39g39syn(7).fvalsI])); % 0.7674
tAPVlearn3=mean(timing3([r39g39syn(9).fvalsA r39g39syn(9).fvalsB r39g39syn(9).fvalsC r39g39syn(9).fvalsI])); % 0.7674
pACbase1=mean(r39g39syn(1).LCtargets);
pACbase2=mean(r39g39syn(3).LCtargets);
tACbase1=mean(timing3([r39g39syn(1).fvalsA r39g39syn(1).fvalsB r39g39syn(1).fvalsC r39g39syn(1).fvalsI])); % 0.7674
tACbase2=mean(timing3([r39g39syn(3).fvalsA r39g39syn(3).fvalsB r39g39syn(3).fvalsC r39g39syn(3).fvalsI])); % 0.7674
%LClearnAC=[r39g39syn(4).LCtargets;r39g39syn(6).LCtargets;r39g39syn(8).LCtargets;r39g39syn(10).LCtargets];
LClearnAC=[r39g39syn(4).LClearn;r39g39syn(6).LClearn;r39g39syn(8).LClearn;r39g39syn(10).LClearn];
   

timesA=[];timesB=[];timesC=[];timesI=[];listchars=[];timechars=[];
    for i=[4 6 8 10] % learning ACSF
        timesA=[timesA timing3(r39g39syn(i).fvalsA)];
        timesB=[timesB timing3(r39g39syn(i).fvalsB)];
        timesC=[timesC timing3(r39g39syn(i).fvalsC)];
        timesI=[timesI timing3(r39g39syn(i).fvalsI)];
    end
    LClearnA1=LClearnAC;
    countA=0;countB=0;countC=0;countI=0;
    listchars=LClearnA1;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
                else if listchars(i)=='i';countI=countI+1;timechars(i)=timesI(countI);
                    end;end;end;end
    end
    TClearnA1=timechars;
    
 runanalysis(1-(LClearnA1=='a' | LClearnA1=='c'),1,1-pACbase,0.005,0)  
figure;hold on;
% Experiment 1 - hit A, starting at baseline
    plot(TClearnA1,1-p05(1:end-1),'.')
    plot(TClearnA1,1-p95(1:end-1),'.')
    plot(TClearnA1,1-pmid(1:end-1))
plot(tACbase1,pACbase1,'b.','Markersize',15)
plot(tACbase2,pACbase2,'b.','Markersize',15)
plot(tAPVbase,pAPVbase,'r.','Markersize',15)
plot(tAPVlearn1,pAPVlearn1,'r.','Markersize',15)
plot(tAPVlearn2,pAPVlearn2,'r.','Markersize',15)
plot(tAPVlearn3,pAPVlearn3,'r.','Markersize',15)













% r37g7
load /cardinal4/r37g7/LearningCurves.mat
%%%% 
% go to     /home/jcharles/matlab/learninganaly/Learning Analysis/IndividualAnalysisEM
% runanalysis(1-(LClearnA1=='a'),1,1-LCbaseA,0.005,0)
% click on  'resultsindividual.mat'
figure;hold on;
% Experiment 1 - hit A, starting at baseline
    subplot(2,2,1);hold on;
    plot(TClearnA1,e1.p05(1:end-1),'.')
    plot(TClearnA1,e1.p95(1:end-1),'.')
    plot(TClearnA1,e1.pmid(1:end-1))
    plot([min(TClearnA1)-20 max(TClearnA1)+20],[1-LCbaseA 1-LCbaseA],'k')
    plot(mean(TClearn_apvA1),mean(1-sum(LClearn_apvA1=='a')/length(LClearn_apvA1)),'+','Markersize',15,'Color','r')
    plot(min(TClearnA1)-5,1-LCapvA,'+','Markersize',15,'Color','r')
    ylim([0 1])
    xlim([2495 2555])
% Experiment 2 - hit B, starting at non-baseline
% runanalysis(1-(LClearnB1=='b'),1,1-exp2preB,0.005,0)
subplot(2,2,2);hold on;
    plot(TClearnB1,e2.p05(1:end-1),'.')
    plot(TClearnB1,e2.p95(1:end-1),'.')
    plot(TClearnB1,e2.pmid(1:end-1))
    plot([min(TClearnB1)-20 max(TClearnB1)+20],[1-LCbaseB 1-LCbaseB],'k') % baseline B
    plot([min(TClearnB1)-20 max(TClearnB1)+20],[1-exp2preB 1-exp2preB],'k') % p(B) before exp 2 (diff from baseline b/c WN on A)
    plot(mean(TClearn_apvB1(12).data),mean(1-sum(LClearn_apvB1(12).data=='b')/length(LClearn_apvB1(12).data)),'+','Markersize',15,'Color','r')
    plot(mean(TClearn_apvB1(14).data),mean(1-sum(LClearn_apvB1(14).data=='b')/length(LClearn_apvB1(14).data)),'+','Markersize',15,'Color','r')
    plot(mean(TClearn_apvB1(16).data),mean(1-sum(LClearn_apvB1(16).data=='b')/length(LClearn_apvB1(16).data)),'+','Markersize',15,'Color','r') 
    plot(min(TClearnB1)-5,1-LCapvB,'+','Markersize',15,'Color','g') % baseline apvB
    plot(min(TClearnB1)-5,mean(1-sum(LClearn_apvA1=='b')/length(LClearn_apvA1)),'+','Markersize',15,'Color','r') % p(B) with apv before exp 2
    xlim([2545 2615])
    ylim([0 1])
% Experiment 3 - hit A, starting near baseline
%runanalysis(1-(LClearnA2=='a'),1,1-exp3preA,0.005,0)
subplot(2,2,3);hold on;
    plot(TClearnA2,e3.p05(1:end-1),'.')
    plot(TClearnA2,e3.p95(1:end-1),'.')
    plot(TClearnA2,e3.pmid(1:end-1))
    plot([min(TClearnA2)-20 max(TClearnA2)+20],[1-LCbaseA 1-LCbaseA],'k') % baseline B
    plot([min(TClearnA2)-20 max(TClearnA2)+20],[1-exp3preA 1-exp3preA],'k') % p(B) before exp 2 (diff from baseline b/c WN on A)
    plot(mean(TClearn_apvA2(18).data),mean(1-sum(LClearn_apvA2(18).data=='a')/length(LClearn_apvA2(18).data)),'+','Markersize',15,'Color','r')
    plot(mean(TClearn_apvA2(20).data),mean(1-sum(LClearn_apvA2(20).data=='a')/length(LClearn_apvA2(20).data)),'+','Markersize',15,'Color','r')
    plot(min(TClearnA2)-5,1-LCapvA,'+','Markersize',15,'Color','g') % baseline apvB
    plot(min(TClearnA2)-5,1-exp3preAapv,'+','Markersize',15,'Color','r') % p(B) with apv before exp 2
    xlim([2615 2700])
    ylim([0 1])
% Experiment 4 - hit B (and then hit A) with AP5 on, starting at non-baseline
%runanalysis(1-(LClearnB2=='b'),1,1-exp4preBapv,0.005,0)
subplot(2,2,4);hold on;
    plot(TClearnB2,e4.p05(1:end-1),'.','Color','r')
    plot(TClearnB2,e4.p95(1:end-1),'.','Color','r')
    plot(TClearnB2,e4.pmid(1:end-1),'Color','r')
    plot([min(TClearnA2)-20 max(TClearnA2)+20],[1-LCbaseB 1-LCbaseB],'k') % baseline B
    plot([min(TClearnA2)-20 max(TClearnA2)+20],[1-exp4preBapv 1-exp4preBapv],'k') % p(B) before exp 2 (diff from baseline b/c WN on A)
    xlim([2672 2680])
    ylim([0 1])

%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
clear all
load /cardinal4/r37g7/SyntaxFF.mat
% Control - huge variability reduction
    % 416 - first inactivation (A targeted so measure B) - 52% reduction
        mean(std(pitch416amB(300:500,:)')) % 37.11
        mean(std(pitch416apvB(300:500,:)')) % 17.89
    % 421 - final inactivation (A targeted so measure B) - 40% reduction
        mean(std(pitch421amB(300:500,:)')) % 35.54
        mean(std(pitch421apvB(300:500,:)')) % 21.55
%


% Control - FF variability reduction and pitch shift reversion
% unclear if there really is learning or var reduction...
figure;hold on;
plot(tv424apvC(1:3:76),pitch424apvC(220,1:3:76),'*','Color','r')
plot(timing3(fvals425C(1:3:160)),pitch425C(220,1:3:160),'*','Color','k')
plot(tv426Catch,pitch426Catch(220,:),'*','Color','b')
plot(tv426apvC,pitch426apvC(220,:),'*','Color','r')



% Step 1: Get times of notes and order of notes for each baseline, AP5, and
% learning folder.
% First experiment hits A
% Order of notes
i=5;
    Experiment1(i).listchars=notestats;
% Run fvals
    Experiment1(i).fvA=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvB=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvD=findwnoteJC('batchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvE=findwnoteJC('batchnotes','e','','',0,[2000 2700],8500,1,'obs0',1);
%%%
    Experiment1(i).listcharscatch=notestats;
    Experiment1(i).fvAcatch=findwnoteJC('batchcatchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvBcatch=findwnoteJC('batchcatchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvDcatch=findwnoteJC('batchcatchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvEcatch=findwnoteJC('batchcatchnotes','e','','',0,[2000 2700],8500,1,'obs0',1);


% Learning experiment 1 - A vs B,D,E
% minor problem - look at March 28 data to replace 5 and 6
                timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];
                    i=1;
                    while Experiment1(i).baseline==1 && i<5
                        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
                        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
                        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
                        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
                        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
                        i=i+1;
                    end
                    LCbase=listchars;
                    countA=0;countB=0;countD=0;countE=0;
                    for i=1:length(listchars)
                        if listchars(i)=='a'
                            countA=countA+1;
                            timechars(i)=timesA(countA);
                        else if listchars(i)=='b'
                            countB=countB+1;
                            timechars(i)=timesB(countB);
                        else if listchars(i)=='d'
                            countD=countD+1;
                            timechars(i)=timesD(countD);
                        else if listchars(i)=='e'
                            countE=countE+1;
                            timechars(i)=timesE(countE);

                        end;end;end;end
                    end
                    TCbase=timechars;
                    windowlength=30;
                    for i=1:floor(length(listchars)/windowlength)
                        window=[(i-1)*windowlength+1:i*windowlength];
                        probNote(i)=sum(listchars(window)=='a')/windowlength;
                        thetime(i)=median(timechars(window));
                    end
                    figure;hold on;plot(thetime,1-probNote,'*','Color','k')
                % pre
                   for i=5:7
                        timesA=[timing3(Experiment1(i).fvAcatch)];
                        timesB=[timing3(Experiment1(i).fvBcatch)];
                        timesD=[timing3(Experiment1(i).fvDcatch)];
                        timesE=[timing3(Experiment1(i).fvEcatch)];
                        listchars=[Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
                        countA=0;countB=0;countD=0;countE=0;
                        timechars=[];
                    for j=1:length(listchars)
                        if listchars(j)=='a'
                            countA=countA+1;
                            timechars(j)=timesA(countA);
                        else if listchars(j)=='b'
                            countB=countB+1;
                            timechars(j)=timesB(countB);
                        else if listchars(j)=='d'
                            countD=countD+1;
                            timechars(j)=timesD(countD);
                        else if listchars(j)=='e'
                            countE=countE+1;
                            timechars(j)=timesE(countE);
                        end;end;end;end
                    end
                        probNote=[];
                        thetime=[];
                        windowlength=30;
                        for j=1:floor(length(listchars)/windowlength)
                            window=[(j-1)*windowlength+1:j*windowlength];
                            probNote(j)=sum(listchars(window)=='a')/windowlength;
                            thetime(j)=median(timechars(window));
                        end
                        if Experiment1(i).acsf==1;
                            plot(mean(thetime)+400,mean(1-probNote),'+','Markersize',15,'Color','b')                            
                        else
                            plot(mean(thetime)+400,mean(1-probNote),'+','Markersize',15,'Color','r')                            
                        end
                end
                
                % learning
                for i=8:22
                        timesA=[timing3(Experiment1(i).fvAcatch)];
                        timesB=[timing3(Experiment1(i).fvBcatch)];
                        timesD=[timing3(Experiment1(i).fvDcatch)];
                        timesE=[timing3(Experiment1(i).fvEcatch)];
                        listchars=[Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
                        countA=0;countB=0;countD=0;countE=0;
                        timechars=[];
                    for j=1:length(listchars)
                        if listchars(j)=='a'
                            countA=countA+1;
                            timechars(j)=timesA(countA);
                        else if listchars(j)=='b'
                            countB=countB+1;
                            timechars(j)=timesB(countB);
                        else if listchars(j)=='d'
                            countD=countD+1;
                            timechars(j)=timesD(countD);
                        else if listchars(j)=='e'
                            countE=countE+1;
                            timechars(j)=timesE(countE);
                        end;end;end;end
                    end
                        probNote=[];
                        thetime=[];
                        windowlength=10;
                        for j=1:floor(length(listchars)/windowlength)
                            window=[(j-1)*windowlength+1:j*windowlength];
                            probNote(j)=sum(listchars(window)=='a')/windowlength;
                            thetime(j)=median(timechars(window));
                        end
                        if Experiment1(i).acsf==1;
                            plot(thetime,1-probNote,'*','Color','b')                            
                        else
                                plot(mean(timechars),(1-sum(listchars=='a')/length(listchars)),'+','Markersize',10,'Color','r')
                        end
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Learning curves %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline ACSF - probability
    listchars=[];
        for i=1:4
            listchars=[listchars;Experiment1(i).listchars(find(Experiment1(i).listchars=='a' | Experiment1(i).listchars=='b' | Experiment1(i).listchars=='d' | Experiment1(i).listchars=='e'))];
        end
        LCbaseA=sum(listchars=='a')/length(listchars); % 0.424
        LCbaseB=sum(listchars=='b')/length(listchars); % 0.286
% Baseline AP5 - probability
    i=7;
    listchars=[];
    listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
        LCapvA=sum(listchars=='a')/length(listchars); % 0.349
        LCapvB=sum(listchars=='b')/length(listchars); % 0.248
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 1 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning ACSF - first experiment - hit A
    timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
    for i=[8 10] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
    end
    LClearnA1=listchars;
    countA=0;countB=0;countD=0;countE=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);
                else if listchars(i)=='e';countE=countE+1;timechars(i)=timesE(countE);
                    end;end;end;end
    end
    TClearnA1=timechars;
% Learning APV - first experiment - hit A
    timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
    for i=9 % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
    end
    LClearn_apvA1=listchars;
    countA=0;countB=0;countD=0;countE=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);
                else if listchars(i)=='e';countE=countE+1;timechars(i)=timesE(countE);
                    end;end;end;end
    end
    TClearn_apvA1=timechars;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 2 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning ACSF - second experiment - hit B

    % preB - since not starting at baseline
    exp2preB=0.4728; % (avg for experiments 8,10)

%%% Learning
    timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
    for i=[11 13 15] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
    end
    LClearnB1=listchars;
    countA=0;countB=0;countD=0;countE=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);
                else if listchars(i)=='e';countE=countE+1;timechars(i)=timesE(countE);
                    end;end;end;end
    end
    TClearnB1=timechars;
% Learning APV - second experiment - hit B
    for i=[12 14 16] % learning ACSF
        timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        LClearn_apvB1(i).data=Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'));
        listchars=LClearn_apvB1(i).data;
        countA=0;countB=0;countD=0;countE=0;
        for j=1:length(listchars)
            if listchars(j)=='a';countA=countA+1;timechars(j)=timesA(countA);
            else if listchars(j)=='b';countB=countB+1;timechars(j)=timesB(countB);
                else if listchars(j)=='d';countD=countD+1;timechars(j)=timesD(countD);
                    else if listchars(j)=='e';countE=countE+1;timechars(j)=timesE(countE);
                        end;end;end;end
        end
        TClearn_apvB1(i).data=timechars;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 3 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Learning pre
        exp3preA=0.3434; % from experiment 15
    % apv pre
        exp3preAapv=mean(sum(LClearn_apvB1(16).data=='a')/length(LClearn_apvB1(16).data)); 
  %%% Learning
    timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
    for i=[17 19 21] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
    end
    LClearnA2=listchars;
    countA=0;countB=0;countD=0;countE=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);
                else if listchars(i)=='e';countE=countE+1;timechars(i)=timesE(countE);
                    end;end;end;end
    end
    TClearnA2=timechars;
% Learning APV - second experiment - hit B
    for i=[18 20] % learning ACSF
        timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        LClearn_apvA2(i).data=Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'));
        listchars=LClearn_apvA2(i).data;
        countA=0;countB=0;countD=0;countE=0;
        for j=1:length(listchars)
            if listchars(j)=='a';countA=countA+1;timechars(j)=timesA(countA);
            else if listchars(j)=='b';countB=countB+1;timechars(j)=timesB(countB);
                else if listchars(j)=='d';countD=countD+1;timechars(j)=timesD(countD);
                    else if listchars(j)=='e';countE=countE+1;timechars(j)=timesE(countE);
                        end;end;end;end
        end
        TClearn_apvA2(i).data=timechars;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 4 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  
% not much data to calculate pre-probabilities - so use setting #2?
exp4preBapv=0.375;

% Learning APV - fourth experiment - hit B -- not enough data from the reverse part
    timesA=[];timesB=[];timesD=[];timesE=[];listchars=[];timechars=[];
    for i=[22] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        timesE=[timesE timing3(Experiment1(i).fvEcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='d' | Experiment1(i).listcharscatch=='e'))];
    end
    LClearnB2=listchars;
    countA=0;countB=0;countD=0;countE=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);
                else if listchars(i)=='e';countE=countE+1;timechars(i)=timesE(countE);
                    end;end;end;end
    end
    TClearnB2=timechars;





%%% binom - reasonable approximation
i=3
listchars=[listchars;Experiment1(i).listchars(find(Experiment1(i).listchars=='a' | Experiment1(i).listchars=='b' | Experiment1(i).listchars=='d' | Experiment1(i).listchars=='e'))];
mchar=(listchars=='a');
pA=sum(mchar)/length(mchar); % 0.33 for i=2, 0.36 for i=3
p2=[];
for i=1:length(mchar)-1
    p2(i)=sum(mchar(i:i+1));
end
sum(p2==0)/length(p2); % 0.40 vs 0.43, 0.36 vs 0.40
sum(p2==1)/length(p2); % 0.52 vs 0.47, 0.56 vs 0.47
sum(p2==2)/length(p2); % 0.07 vs 0.10, 0.08 vs 0.13
% These are as would be predicted.
mchar=(listchars=='a');
pA=sum(mchar)/length(mchar); % 0.47 for i=3, 0.4 for i=5
p2=[];
for i=1:length(mchar)-1
    p2(i)=sum(mchar(i:i+1));
end
sum(p2==0)/length(p2); % 0.29 vs 0.28   % 0.34 for i=5 vs. 0.36
sum(p2==1)/length(p2); % 0.49 vs 0.5   % 0.51 for i=5 vs. 0.48
sum(p2==2)/length(p2); % 0.22 vs 0.22  % 0.15 for i=5 vs. 0.16
% These are as would be predicted.
mchar=(listchars=='c');
pC=sum(mchar)/length(mchar); % 0.47 for i=3, 0.4 for i=5
p2=[];
for i=1:length(mchar)-1
    p2(i)=sum(mchar(i:i+1));
end
sum(p2==0)/length(p2); % 0.29 vs 0.28   % 0.34 for i=5 vs. 0.36
sum(p2==1)/length(p2); % 0.49 vs 0.5   % 0.51 for i=5 vs. 0.48
sum(p2==2)/length(p2); % 0.22 vs 0.22  % 0.15 for i=5 vs. 0.16
% These are as would be predicted.

%%%%%%%%
%%%%%%%%
%%%%%%%%


% pu3bk24
load /cardinal6/pu3bk24/LearningCurves.mat
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
        % baseline - probability of transition at baseline
        cd /home/jcharles/matlab/learninganaly/Learning Analysis/IndividualAnalysisEM
        runanalysis(1-(LClearn=='d'),1,1-pBaseD,0.005,0)
        % 0.005 is the default setting
        % 0 means that the state starts at the baseline probability instead of
            % using the first few learning trials to estimate the current state.
            % This makes sense with our data because we have much more baseline
            % data than learning data.  It would not make sense if we had no
            % baseline data and were just assuming a priori that p=0.5 (chance).

        % click on  'resultsindividual.mat'
        figure;hold on;
        % Experiment 1 - hit D, starting at baseline
            plot(TClearn,p05(1:end-1),'.')
            plot(TClearn,p95(1:end-1),'.')
            plot(TClearn,pmid(1:end-1))
            plot([min(TClearn)-20 max(TClearn)+20],[1-pBaseD 1-pBaseD],'k')
            plot(mean(TClearn_apv),mean(1-sum(LClearn_apv=='d')/length(LClearn_apv)),'+','Markersize',15,'Color','r')
            plot(min(TClearn)-5,1-pBaseDapv,'+','Markersize',15,'Color','r')
            ylim([0 1])
            xlim([8255 8325])

        %%% CV controls
        mean(std(pitch1213acsfCallforFF(220:280,:)'))/mean(mean(pitch1213acsfCallforFF(220:280,:)'))
        mean(std(pitch1213apvCallforFF(220:280,:)'))/mean(mean(pitch1213apvCallforFF(220:280,:)'))
        mean(std(pitchpreC(220:280,:)'))/mean(mean(pitchpreC(220:280,:)'))
        mean(std(pitch1210apvC(220:280,:)'))/mean(mean(pitch1210apvC(220:280,:)'))


%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%

% Baseline
    pBaseD=0.5578;
    pBaseDapv=0.5441;
% Learning ACSF
    timesC=[];timesD=[];listchars=[];
    for i=[6 7] % learning ACSF
        timesC=[timesC timing3(Experiment1(i).fvCcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='c' | Experiment1(i).listcharscatch=='d'))];
    end
    LClearn=listchars;
    countC=0;countD=0;
    timechars=[];
    for i=1:length(listchars)
        if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
        else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);;
        end;end;
    end
    TClearn=timechars;
% Learning APV
        i=8;
        timesC=[];timesD=[];
        timesC=[timesC timing3(Experiment1(i).fvCcatch)];
        timesD=[timesD timing3(Experiment1(i).fvDcatch)];
        listchars=[Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='c' | Experiment1(i).listcharscatch=='d'))];
    LClearn_apv=listchars;
    countC=0;countD=0;
    timechars=[];
    for i=1:length(listchars)
        if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
        else if listchars(i)=='d';countD=countD+1;timechars(i)=timesD(countD);;
        end;end;
    end
    TClearn_apv=timechars;




 % Step 1: Get times of notes and order of notes for each baseline, AP5, and
% learning folder.
% First experiment hits A
% Order of notes
i=3;
   Experiment1(i).listchars=notestats;

    Experiment1(i).fvC=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvD=findwnoteJC('batchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
%%%
    Experiment1(i).listcharscatch=notestats;
    Experiment1(i).fvCcatch=findwnoteJC('batchcatchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvDcatch=findwnoteJC('batchcatchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);

%
figure;hold on;
for i=3:5
    timesC=[timing3(Experiment1(i).fvC)];
    timesD=[timing3(Experiment1(i).fvD)];
    listchars=[Experiment1(i).listchars(find(Experiment1(i).listchars=='c' | Experiment1(i).listchars=='d'))];
    countC=0;countD=0;
    timechars=[];
    for j=1:length(listchars)
        if listchars(j)=='c'
            countC=countC+1;
            timechars(j)=timesC(countC);
        else if listchars(j)=='d'
                countD=countD+1;
                timechars(j)=timesD(countD);
            end;end
    end
    probNote=[];
    thetime=[];
    if Experiment1(i).acsf==1
        windowlength=50;
        for j=1:floor(length(listchars)/windowlength)
            window=[(j-1)*windowlength+1:j*windowlength];
            probNote(j)=sum(listchars(window)=='d')/windowlength;
            thetime(j)=median(timechars(window));
        end
        plot(thetime,(1-probNote),'*','Color','b')
    else
        windowlength=50;
        for j=1:floor(length(listchars)/windowlength)
            window=[(j-1)*windowlength+1:j*windowlength];
            probNote(j)=sum(listchars(window)=='d')/windowlength;
            thetime(j)=median(timechars(window));
        end

        plot((thetime),(1-probNote),'*','Color','r')
    end
end

% learning
for i=6:8
    timesC=[timing3(Experiment1(i).fvCcatch)];
    timesD=[timing3(Experiment1(i).fvDcatch)];
    listchars=[Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='c' | Experiment1(i).listcharscatch=='d'))];
    countC=0;countD=0;
    timechars=[];
    for j=1:length(listchars)
        if listchars(j)=='c'
            countC=countC+1;
            timechars(j)=timesC(countC);
        else if listchars(j)=='d'
                countD=countD+1;
                timechars(j)=timesD(countD);
            end;end
    end
    probNote=[];
    thetime=[];
    windowlength=10;
    for j=1:floor(length(listchars)/windowlength)
        window=[(j-1)*windowlength+1:j*windowlength];
        probNote(j)=sum(listchars(window)=='d')/windowlength;
        thetime(j)=median(timechars(window));
    end
    if Experiment1(i).acsf==1
        plot(thetime,(1-probNote),'*','Color','k')
    else
        plot(mean(thetime),mean(1-probNote),'+','Markersize',15,'Color','r')
    end
end

%%% binomial assumption - reasonable approximation
i=5;
listchars=[Experiment1(i).listchars(find(Experiment1(i).listchars=='c' | Experiment1(i).listchars=='d'))];
mchar=(listchars=='c');
pC=sum(mchar)/length(mchar); % 0.47 for i=3, 0.4 for i=5
p2=[];
for i=1:length(mchar)-1
    p2(i)=sum(mchar(i:i+1));
end
sum(p2==0)/length(p2); % 0.29 vs 0.28   % 0.34 for i=5 vs. 0.36
sum(p2==1)/length(p2); % 0.49 vs 0.5   % 0.51 for i=5 vs. 0.48
sum(p2==2)/length(p2); % 0.22 vs 0.22  % 0.15 for i=5 vs. 0.16
% These are as would be predicted.



%%%%%%%%
%%%%%%%%
%%%%%%%%

% pu87bk30
load /cardinal6/pu87bk30/LearningCurves.mat
        cd /home/jcharles/matlab/learninganaly/Learning Analysis/IndividualAnalysisEM

        %runanalysis(1-(LClearn=='a'),1,1-pBaseA,0.005,0)
        % click on  'resultsindividual.mat'
% Learning
figure;hold on;
        % Experiment 1 - hit D, starting at baseline
            plot(TClearn,p05(1:end-1),'.')
            plot(TClearn,p95(1:end-1),'.')
            plot(TClearn,pmid(1:end-1))
            plot([240 320],[1-pBaseA 1-pBaseA],'k')
            plot(mean(TClearnAPV),mean(1-sum(LClearnAPV=='a')/length(LClearnAPV)),'+','Markersize',15,'Color','r')
            plot(min(TClearn)-5,1-pBaseAapv,'+','Markersize',15,'Color','r')
            ylim([0 1])
            xlim([245 315])
% Branch point
figure;imagesc([0:0.0032:0.0032*554],f,log(avZ));syn;ylim([0,1e4]);xlim([0.3 1.7])
% Variability reduction as metric of effectiveness - GREAT0.
    mean(std(pitch113apvA(200:300,:)'))/mean(mean(pitch113apvA(200:300,:)'))
    mean(std(pitch109apvA(200:300,:)'))/mean(mean(pitch109apvA(200:300,:)'))
    mean(std(pitchPRE(200:300,:)'))/mean(mean(pitchPRE(200:300,:)'))
    mean(std(pitch111acsfA(200:300,:)'))/mean(mean(pitch111acsfA(200:300,:)'))





% Baseline
pBaseA=0.7565;
pBaseAapv=0.7671;

% Learning ACSF
    timesA=[];timesB=[];listchars=[];
    for i=[4 5] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
    end
    LClearn=listchars;
    countA=0;countB=0;
    timechars=[];
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
        end;end;
    end
    TClearn=timechars;
% Learning APV
    timesA=[];timesB=[];listchars=[];
    for i=6 % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
    end
    LClearnAPV=listchars;
    countA=0;countB=0;
    timechars=[];
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
        end;end;
    end
    TClearnAPV=timechars;
%%%%



 % Step 1: Get times of notes and order of notes for each baseline, AP5, and
% learning folder.
% First experiment hits A
% Order of notes
i=7;
   Experiment3(i).listchars=notestats;

    Experiment3(i).fvA=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment3(i).fvB=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
%%%
    Experiment3(i).listcharscatch=notestats;
    Experiment3(i).fvAcatch=findwnoteJC('batchcatchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment3(i).fvBcatch=findwnoteJC('batchcatchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    figure;hold on
for i=3
    timesA=[timing3(Experiment1(i).fvAcatch)];
    timesB=[timing3(Experiment1(i).fvBcatch)];
    listchars=[Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
    countA=0;countB=0;
    timechars=[];
    for j=1:length(listchars)
        if listchars(j)=='a'
            countA=countA+1;
            timechars(j)=timesA(countA);
        else if listchars(j)=='b'
                countB=countB+1;
                timechars(j)=timesB(countB);
            end;end
    end
    probNote=[];
    thetime=[];
    windowlength=10;
    for j=1:floor(length(listchars)/windowlength)
        window=[(j-1)*windowlength+1:j*windowlength];
        probNote(j)=sum(listchars(window)=='a')/windowlength;
        thetime(j)=median(timechars(window));
    end
    if Experiment1(i).acsf==1
        plot(thetime,(1-probNote),'*','Color','k')
    else
        plot(mean(thetime),mean(1-probNote),'+','Markersize',15,'Color','r')
    end
end

i=5;
listchars=[Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
mchar2=(listchars=='b');


% bk75bk62

        %runanalysis(1-(LClearn=='a'),1,1-pBaseA,0.005,0)
        % click on  'resultsindividual.mat'

i=1;
   Experiment3(i).listchars=notestats;

    Experiment3(i).fvB=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment3(i).fvC=findwnoteJC('batchnotes','a','cc','',0,[2000 2700],8500,1,'obs0',1);
%%%
i=4
    Experiment3(i).listcharscatch=notestats;
    Experiment3(i).fvBcatch=findwnoteJC('batchcatchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment3(i).fvCcatch=findwnoteJC('batchcatchnotes','a','c-','',0,[2000 2700],8500,1,'obs0',1);


% Baseline
pBaseB=0.6027;
pBaseBapv=0.6538; % from 927_2mMapv
% Experiment 1- Learning ACSF
    timesB=[];timesC=[];listchars=[];
    for i=2 % learning ACSF
        timesB=[timesB timing4(Experiment1(i).fvBcatch)];
        timesC=[timesC timing4(Experiment1(i).fvCcatch)];
        for j=3:length(Experiment1(i).listcharscatch)
            if Experiment1(i).listcharscatch(j)=='a' & Experiment1(i).listcharscatch(j-1)=='b'
                listchars=[listchars;Experiment1(i).listcharscatch(j-1)];
            else if Experiment1(i).listcharscatch(j)=='a' & Experiment1(i).listcharscatch(j-2)=='c'
                    listchars=[listchars;Experiment1(i).listcharscatch(j-2)];
                end
            end
        end
    end
    LClearn=listchars;
    countB=0;countC=0;
    TClearn=[];
    timechars=[];
    for i=1:length(listchars)
        if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
        else if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
        end;end;
    end
    TClearn=timechars;
% Learning APV
    timesB=[];timesC=[];listchars=[];
    for i=4 % learning APV
        timesB=[timesB timing4(Experiment1(i).fvBcatch)];
        timesC=[timesC timing4(Experiment1(i).fvCcatch)];
        for j=3:length(Experiment1(i).listcharscatch)
            if Experiment1(i).listcharscatch(j)=='a' & Experiment1(i).listcharscatch(j-1)=='b'
                listchars=[listchars;Experiment1(i).listcharscatch(j-1)];
            else if Experiment1(i).listcharscatch(j)=='a' & Experiment1(i).listcharscatch(j-2)=='c'
                    listchars=[listchars;Experiment1(i).listcharscatch(j-2)];
                end
            end
        end
    end
    LClearnAPV=listchars;
    countB=0;countC=0;
    timechars=[];
    for i=1:length(listchars)
        if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
        else if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
        end;end;
    end
    TClearnAPV=timechars;
%%%%
        /home/jcharles/matlab/learninganaly/Learning Analysis/IndividualAnalysisEM
        runanalysis(1-(LClearn=='b'),1,1-pBaseB,0.005,0)
        runanalysis(1-(LClearnAPV=='b'),1,1-pBaseBapv,0.005,0)
% Learning ACSF
figure;hold on;
        % Experiment 1 - hit BA, starting at baseline, ACSF on
            plot(TClearn,e1.p05(1:end-1),'.')
            plot(TClearn,e1.p95(1:end-1),'.')
            plot(TClearn,e1.pmid(1:end-1))
            plot([7040 7080],[1-pBaseB 1-pBaseB],'k')
            ylim([0 1])
        % Experiment 2 - hit BA, starting at baseline, APV on
        % APV clearly blocks much learning, but ACSF+WN at beginning for
        % 30min
            plot(TClearnAPV-420,e2.p05(1:end-1),'.','Color','r')
            plot(TClearnAPV-420,e2.p95(1:end-1),'.','Color','r')
            plot(TClearnAPV-420,e2.pmid(1:end-1),'Color','r')
            plot([7040 7080],[1-pBaseBapv 1-pBaseBapv],'r')
            ylim([0 1])
                               xlim([7040 7060])
            
            
% bk13bk12
% C is low stack
i=3;
   Experiment1(i).listchars=notestats;

    Experiment1(i).fvA=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvB=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
%%%
i=8
    Experiment1(i).listcharscatch=notestats;
    Experiment1(i).fvAcatch=findwnoteJC('batchcatchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvBcatch=findwnoteJC('batchcatchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);

% Learning
    pBaseB=0.4040;
    pBaseBapv=0.3839;
    pDuringBapv=0.3636;
    pBaseA=0.5960;
    pBaseA2=0.2059; % p before Exp2
    pBaseAapv=0.6161;
    pBaseA2apv=0.6364; % p before Exp2
  % Experiment 1- Learning hitting B
        timesA=[];timesB=[];listchars=[];
        for i=[3 5] % learning ACSF
            timesA=[timesA timing4(Experiment1(i).fvAcatch)];
            timesB=[timesB timing4(Experiment1(i).fvBcatch)];
            listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
        end
        LClearnB=listchars;
        countA=0;countB=0;
        for i=1:length(listchars)
            if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
            else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            end;end;
        end
        TClearnA=timechars;
        % Reversion APV
            timesA=[];timesB=[];listchars=[];
            for i=4 % learning ACSF
                timesA=[timesA timing4(Experiment1(i).fvAcatch)];
                timesB=[timesB timing4(Experiment1(i).fvBcatch)];
                listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
            end
            LClearnAPV=listchars;
            countA=0;countB=0;
            timechars=[];
            for i=1:length(listchars)
                if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
                else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
                end;end;
            end
            TClearnAPV=timechars;
  % Experiment 2- Learning hitting A
        timesA=[];timesB=[];listchars=[];
        for i=[6] % learning ACSF
            timesA=[timesA timing4(Experiment1(i).fvAcatch)];
            timesB=[timesB timing4(Experiment1(i).fvBcatch)];
            listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
        end
        LClearnB=listchars;
        countA=0;countB=0;
        for i=1:length(listchars)
            if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
            else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            end;end;
        end
        TClearnB=timechars;
        % Reversion APV
            timesA=[];timesB=[];listchars=[];
            for i=7 % learning ACSF
                timesA=[timesA timing4(Experiment1(i).fvAcatch)];
                timesB=[timesB timing4(Experiment1(i).fvBcatch)];
                listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b'))];
            end
            LClearnAPV=listchars;
            countA=0;countB=0;
            timechars=[];
            for i=1:length(listchars)
                if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
                else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
                end;end;
            end
            TClearnAPV2=timechars;

        cd /home/jcharles/matlab/learninganaly/LearningAnalysis/IndividualAnalysisEM
       % runanalysis(1-(LClearnA=='b'),1,1-pBaseB,0.005,0)
      % runanalysis((LClearnB=='a'),1,pBaseA2,0.005,0)
figure;hold on;
        % Experiment 1 - hit B, starting at baseline, ACSF on
            plot(TClearnA,e1.p05(1:end-1),'.')
            plot(TClearnA,e1.p95(1:end-1),'.')
            plot(TClearnA,e1.pmid(1:end-1))
            plot([8300 8370],[1-pBaseB 1-pBaseB],'k')
            plot(mean(TClearnAPV),1-pDuringBapv,'+','Markersize',15,'Color','r')
            plot(min(TClearnA)-5,1-pBaseBapv,'+','Markersize',15,'Color','r')
            ylim([0 1])
        % Experiment 1 - hit A for reversion and beyond
        % Experiment 1 - hit B, starting at baseline, ACSF on
            plot(TClearnB,e2.p05(1:end-1),'.')
            plot(TClearnB,e2.p95(1:end-1),'.')
            plot(TClearnB,e2.pmid(1:end-1))
            plot([8340 8370],[1-pBaseA2 1-pBaseA2],'k')
            plot(mean(TClearnAPV2),0.2692,'+','Markersize',15,'Color','r')
            ylim([0 1])
            xlim([8305 8370])

            
            
            
%%%%%%%%%            
% pu67bk2
%%%%%%%%5
 runanalysis(1-(LClearn=='a'),1,1-pBaseA,0.005,0)
 
         figure;hold on;
        % Experiment 1 - hit A, starting at baseline
            plot(TClearn,e1.p05(1:end-1),'.')
            plot(TClearn,e1.p95(1:end-1),'.')
            plot(TClearn,e1.pmid(1:end-1))
            plot([min(TClearn)-20 max(TClearn)+20],[1-pBaseA 1-pBaseA],'k')
            for i=[4 6 8]
            plot(mean(TClearn_apv(i).data),mean(1-sum(LClearn_apv(i).data=='a')/length(LClearn_apv(i).data)),'+','Markersize',15,'Color','r')
            end
            plot(min(TClearn)-5,1-pBaseAapv,'+','Markersize',15,'Color','r')
            ylim([0 1])
            xlim([4240 4350])            


% SYNSHIFT
% /syn_627_wnon
i=8;
%     Experiment1(i).listchars=notestats;
%         Experiment1(i).fvA=findwnoteJC('batchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
%     %Experiment1(i).fvB=findwnoteJC('batchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
%     Experiment1(i).fvD=findwnoteJC('batchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
%     Experiment1(i).fvE=findwnoteJC('batchnotes','e','','',0,[2000 2700],8500,1,'obs0',1);
%     Experiment1(i).fvC=findwnoteJC('batchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);

% Run fvals
%%%
    Experiment1(i).listcharscatch=notestats;
    Experiment1(i).fvAcatch=findwnoteJC('batchcatchnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvBcatch=findwnoteJC('batchcatchnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvDcatch=findwnoteJC('batchcatchnotes','d','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvEcatch=findwnoteJC('batchcatchnotes','e','','',0,[2000 2700],8500,1,'obs0',1);
    Experiment1(i).fvCcatch=findwnoteJC('batchcatchnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
% Learning pre
        pBaseA=0.952; % from experiment 1
    % apv pre
        pBaseAapv=0.921; % from experiment 2
  %%% Learning
    timesA=[];timesB=[];timesC=[];listchars=[];timechars=[];
    for i=[3 5 7] % learning ACSF
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesC=[timesC timing3(Experiment1(i).fvCcatch)];
        listchars=[listchars;Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='c'))];
    end
    LClearn=listchars;
    countA=0;countB=0;countC=0;
    for i=1:length(listchars)
        if listchars(i)=='a';countA=countA+1;timechars(i)=timesA(countA);
        else if listchars(i)=='b';countB=countB+1;timechars(i)=timesB(countB);
            else if listchars(i)=='c';countC=countC+1;timechars(i)=timesC(countC);
                    end;end;end;
    end
    TClearn=timechars;
% Learning APV - second experiment - hit B
    for i=[4 6 8] % learning ACSF
        timesA=[];timesB=[];timesC=[];listchars=[];timechars=[];
        timesA=[timesA timing3(Experiment1(i).fvAcatch)];
        timesB=[timesB timing3(Experiment1(i).fvBcatch)];
        timesC=[timesC timing3(Experiment1(i).fvCcatch)];
        LClearn_apv(i).data=Experiment1(i).listcharscatch(find(Experiment1(i).listcharscatch=='a' | Experiment1(i).listcharscatch=='b' | Experiment1(i).listcharscatch=='c'));
        listchars=LClearn_apv(i).data;
        countA=0;countB=0;countC=0;
        for j=1:length(listchars)
            if listchars(j)=='a';countA=countA+1;timechars(j)=timesA(countA);
            else if listchars(j)=='b';countB=countB+1;timechars(j)=timesB(countB);
                else if listchars(j)=='c';countC=countC+1;timechars(j)=timesC(countC);
                        end;end;end;
        end
        TClearn_apv(i).data=timechars;
    end

 %%%% 
 runanalysis(1-(LClearn=='a'),1,1-pBaseA,0.005,0)
 
         figure;hold on;
        % Experiment 1 - hit D, starting at baseline
            plot(TClearn,p05(1:end-1),'.')
            plot(TClearn,p95(1:end-1),'.')
            plot(TClearn,pmid(1:end-1))
            plot([min(TClearn)-20 max(TClearn)+20],[1-pBaseA 1-pBaseA],'k')
            for i=[4 6 8]
            plot(mean(TClearn_apv(i).data),mean(1-sum(LClearn_apv(i).data=='a')/length(LClearn_apv(i).data)),'+','Markersize',15,'Color','r')
            end
            plot(min(TClearn)-5,1-pBaseDapv,'+','Markersize',15,'Color','r')
            ylim([0 1])
            xlim([8255 8325])            
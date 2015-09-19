% See Syntax071510.m for calculations
        cd /home/jcharles/matlab/learninganaly/LearningAnalysis/IndividualAnalysisEM
         runanalysis(1-(LClearn=='d'),1,1-pBaseD,0.005,0)
        % 0.005 is the default setting
        % 0 means that the state starts at the baseline probability instead of
            % using the first few learning trials to estimate the current state.
            % This makes sense with our data because we have much more baseline
            % data than learning data.  It would not make sense if we had no
            % baseline data and were just assuming a priori that p=0.5 (chance).

        % click on  'resultsindividual.mat'

        
% summary for plotting
        load /cardinal3/SyntaxLearningFigs/birds.mat
%         % For each inactivation point, information about time, acsf p(transition), apv p(transition), and confidence intervals
%         i=8;
%         bird(i).name='r37g7_exp3';
%         bird(i).days=[0 2];
%         bird(i).acsfmean=[0.6566 0.8686];
%         bird(i).acsf025=[];
%         bird(i).acsf975=[]; % only at baseline - 95% Confidence
%         bird(i).acsf05=[0.6345 0.7910];
%         bird(i).acsf95=[0.6781 0.9202];
%         bird(i).apvmean=[0.7059 0.8750];
%         bird(i).apv05=[0.5524 0.6562];
%         bird(i).apv95=[0.8309 0.9773];
%         bird(i).apv025=[0.5252 0.6165];
%         bird(i).apv975=[0.8490 0.9845];
% 
%         %         [phat,pci]=binofit(sum(LClearnA1=='a'),length(LClearnAPV),0.1)   
%         %         e3.pmid(max(find(TClearnA2-min(TClearn_apvA2)<0)))
figure;subplot(121);hold on;
for i=1:8
    for ii=1:length(bird(i).days)
        if bird(i).days(ii)==0
            % plot(bird(i).days(ii),bird(i).apvmean(ii)-bird(i).acsfmean(ii),'*')
        else
            plot(bird(i).days(ii)+0.2,bird(i).apvmean(ii)-bird(i).apvmean(1),'.','Markersize',15,'Color','r')
            plot(bird(i).days(ii)-0.2,bird(i).acsfmean(ii)-bird(i).acsfmean(1),'.','Markersize',15,'Color','b')
            if i==1
            plot([bird(i).days(ii)-0.2;bird(i).days(ii)+0.2],[bird(i).acsfmean(ii)-bird(i).acsfmean(1);bird(i).apvmean(ii)-bird(i).apvmean(1)],'g')
            else
            plot([bird(i).days(ii)-0.2;bird(i).days(ii)+0.2],[bird(i).acsfmean(ii)-bird(i).acsfmean(1);bird(i).apvmean(ii)-bird(i).apvmean(1)],'k')
            end
        end
    end
end
subplot(122);hold on;
for i=1:8
    for ii=1:length(bird(i).days)
        if bird(i).days(ii)==0
            % plot(bird(i).days(ii),bird(i).apvmean(ii)-bird(i).acsfmean(ii),'*')
        else
            plot(bird(i).days(ii)+0.2,bird(i).apvmean(ii)-bird(i).acsfmean(1),'.','Markersize',15,'Color','r')
            plot(bird(i).days(ii)-0.2,bird(i).acsfmean(ii)-bird(i).acsfmean(1),'.','Markersize',15,'Color','b')
            if i==1
              plot([bird(i).days(ii)-0.2;bird(i).days(ii)+0.2],[bird(i).acsfmean(ii)-bird(i).acsfmean(1);bird(i).apvmean(ii)-bird(i).acsfmean(1)],'g')   
            else
              plot([bird(i).days(ii)-0.2;bird(i).days(ii)+0.2],[bird(i).acsfmean(ii)-bird(i).acsfmean(1);bird(i).apvmean(ii)-bird(i).acsfmean(1)],'k')   
            end

        end
    end
end


 %%
 %%
 %%
        
        
%%%%%%%%%%
%%%%%%%%%%
% pu67bk2
%%%%%%%%%%
%%%%%%%%%%

     load /cardinal4/pu67bk2/LearningCurves.mat
         figure;hold on;
        % Experiment 1 - hit A, starting at baseline (vs sum(A,B,C))
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
mean(std(pitchacsfD(265:320,:)'))
mean(std(pitchapvD(265:320,:)'))

%%%%%%%%%%
%%%%%%%%%%
% pu87bk30
%%%%%%%%%%
%%%%%%%%%%
    load /cardinal6/pu87bk30/LearningCurves.mat
    % Learning
    figure;hold on;
    subplot(211);hold on;
            % Experiment 1 - hit A, starting at baseline
                plot(TClearn,p05(1:end-1),'.')
                plot(TClearn,p95(1:end-1),'.')
                plot(TClearn,pmid(1:end-1))
                plot([240 320],[1-pBaseA 1-pBaseA],'k')
                plot(mean(TClearnAPV),mean(1-sum(LClearnAPV=='a')/length(LClearnAPV)),'+','Markersize',15,'Color','r')
                plot(min(TClearn)-5,1-pBaseAapv,'+','Markersize',15,'Color','r')
                ylim([0 1])
                xlim([245 315])
    % Branch point
        subplot(212);clim=[6 15];imagesc([0:0.0032:0.0032*554],f,log(avZ),clim);syn;ylim([0,1e4]);xlim([0.3 1.7])
    % Variability reduction as metric of effectiveness 
       % 35% reduction --- CVratioEarly=[mean(std(pitch109apvA(200:300,:)'))/mean(mean(pitch109apvA(200:300,:)'))]/[mean(std(pitchPRE(200:300,:)'))/mean(mean(pitchPRE(200:300,:)'))];
       % 53% reduction --- CVratioLate=[mean(std(pitch113apvA(200:300,:)'))/mean(mean(pitch113apvA(200:300,:)'))]/[mean(std(pitch111acsfA(200:300,:)'))/mean(mean(pitch111acsfA(200:300,:)'))];

%%%%%%%%%%
%%%%%%%%%%
% pu3bk24
%%%%%%%%%%
%%%%%%%%%%
    load /cardinal6/pu3bk24/LearningCurves.mat
    % Learning
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
    % Branch point
        figure;imagesc(t,f,log(avZ));syn;ylim([0,1e4]);xlim([-1.3 0.2])
    % Variability reduction as metric of effectiveness
             % 30% reduction ---   CVratioEarly=[mean(std(pitch1213apvCallforFF(220:280,:)'))/mean(mean(pitch1213apvCallforFF(220:280,:)'))]/[mean(std(pitch1213acsfCallforFF(220:280,:)'))/mean(mean(pitch1213acsfCallforFF(220:280,:)'))];
             % 40% reduction ---   CVratioLate=[mean(std(pitch1210apvC(220:280,:)'))/mean(mean(pitch1210apvC(220:280,:)'))]/[mean(std(pitchpreC(220:280,:)'))/mean(mean(pitchpreC(220:280,:)'))];

%%%%%%%%%%
%%%%%%%%%%
% r37g7
%%%%%%%%%%
%%%%%%%%%%
    load /cardinal4/SyntaxAPV/r37g7/LearningCurves.mat  
    load /cardinal4/SyntaxAPV/r37g7/ExperimentFiles.mat
    % Learning
        figure;hold on;
        % Experiment 1 - hit A, starting at baseline   % runanalysis(1-(LClearnA1=='a'),1,1-LCbaseA,0.005,0)

            unidays=find(diff(TClearnA1)>8);
            unistart=[1 unidays+1];
            uniend=[unidays length(TClearnA1)];
            for i=1:length(unistart)
                xvls=[TClearnA1(unistart(i):uniend(i)),TClearnA1(uniend(i):-1:unistart(i))];
                yvls=[1-e1.p05(unistart(i):uniend(i)),1-e1.p95(uniend(i):-1:unistart(i))];
                fill(xvls,yvls,'b')
                plot(TClearnA1(unistart(i):uniend(i)),1-e1.pmid(unistart(i):uniend(i)),'r')
            end
            xv=[mean(timing3(Experiment1(5).fvAcatch)) mean(timing3(Experiment1(7).fvAcatch)) mean(timing3(Experiment1(4).fvAcatch)) 0 mean(timing3(Experiment1(9).fvAcatch))];
            xv=[2470 2480 2490 xv(4:5)];
            for i=[2 5]
                xvls=[xv(i) xv(i)];
                yvls=[probs05{i} probs95{i}];
                plot(xvls,yvls,'r-','Linewidth',2)
                plot(mean(xvls),mean(yvls),'r.','Markersize',15)
            end
            for i=[1 3]
                xvls=[xv(i) xv(i)];
                yvls=[probs05{i} probs95{i}];
                plot(xvls,yvls,'b-','Linewidth',2)
                plot(mean(xvls),mean(yvls),'b.','Markersize',15)
            end
            plot([min(TClearnA1)-50 max(TClearnA1)+20],[mean([probs{1} probs{3}]) mean([probs{1} probs{3}]) ],'k')
            ylim([0 1])
            xlim([2460 2555])
        % Experiment 2 - hit B, starting at non-baseline    % runanalysis(1-(LClearnB1=='b'),1,1-exp2preB,0.005,0)
        subplot(2,3,2);hold on;
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
        % Experiment 3 - hit A, starting near baseline    % runanalysis(1-(LClearnA2=='a'),1,1-exp3preA,0.005,0)
        subplot(2,3,3);hold on;
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
        % Experiment 4 - hit B (and then hit A) with AP5 on, starting at non-baseline    %runanalysis(1-(LClearnB2=='b'),1,1-exp4preBapv,0.005,0)
        subplot(2,3,4);hold on;
            plot(TClearnB2,e4.p05(1:end-1),'.','Color','r')
            plot(TClearnB2,e4.p95(1:end-1),'.','Color','r')
            plot(TClearnB2,e4.pmid(1:end-1),'Color','r')
            plot([min(TClearnA2)-20 max(TClearnA2)+20],[1-LCbaseB 1-LCbaseB],'k') % baseline B
            plot([min(TClearnA2)-20 max(TClearnA2)+20],[1-exp4preBapv 1-exp4preBapv],'k') % p(B) before exp 2 (diff from baseline b/c WN on A)
            xlim([2672 2680])
            ylim([0 1])
    % Branch point
        subplot(2,3,5);
        clim=[6 15];
        figure;imagesc(t,f,log(avZ),clim);syn;ylim([0,1e4]);xlim([-3.4 0])
    % Variability reduction as metric of effectiveness
            % 416 - first inactivation (A targeted so measure B) - 52% reduction
                mean(std(pitch416amB(300:500,:)')) % 37.11
                mean(std(pitch416apvB(300:500,:)')) % 17.89
            % 421 - final inactivation (A targeted so measure B) - 40% reduction
                mean(std(pitch421amB(300:500,:)')) % 35.54
                mean(std(pitch421apvB(300:500,:)')) % 21.55
        % Control - FF variability reduction and pitch shift reversion
        % unclear if there really is learning or var reduction...
                    figure;hold on;
                    plot(tv424apvC(1:3:76),pitch424apvC(220,1:3:76),'*','Color','r')
                    plot(timing3(fvals425C(1:3:160)),pitch425C(220,1:3:160),'*','Color','k')
                    plot(tv426Catch,pitch426Catch(220,:),'*','Color','b')
                    plot(tv426apvC,pitch426apvC(220,:),'*','Color','r')
%%%%%%%%%%
%%%%%%%%%%
% bk75bk62
%%%%%%%%%%
%%%%%%%%%%
    % Looks like APV blocks learning, if one assumes that the earliest one hour
        % is due to ACSF + WN (since JC/TW turned WN on and APV on at the same time
        % instead of waiting for the APV effect).
    % Could this be due to targeting of a convergence point instead of a note after a
        % divergence point.
    % Could this be due the control being before surgery and tethering
        % (note that both have exactly the same number of branch point
        % performances and are approximately the same number of hours - 6.0)
        
     load /cardinal/bk75w62/LearningCurves.mat  
        % Learning ACSF
        figure;subplot(212);hold on;
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
        % Variability reduction - 34%
            mean(std(pitchACSF1001(200:500,:)'))
            mean(std(pitchAPV1001(200:500,:)'))
        % Convergence point
        clim=[4 15];
        subplot(211);imagesc(t,f,log(avZ),clim);syn;ylim([0 1e4]);xlim([-1.2 0.25])
            
            
%%%%%%%%%%
%%%%%%%%%%
% bk13bk12
%%%%%%%%%%
%%%%%%%%%%
load /cardinal/bk13bk12/LearningCurves.mat
% Two experiments in consecutive order - branch point goes to A or B
figure;hold on;subplot(212);hold on;
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
       % Variability reduction
         % baseline (1104tmptestJC - morning b/f WN)
         mean(std(pitchacsfA(200:300,:)'))
         mean(std(pitchacsfB(250:500,:)'))
         % Experiment 1 apv - note A (non-targeted) - 32% reduction
         mean(std(pitchapv1A(200:300,:)'))
         % Experiment 2 apv - note B (non-targeted) - 36% reduction
         mean(std(pitchapv2B(250:500,:)'))
       % Branch point
       clim=[4 15];
        subplot(211);imagesc(t,f,log(avZ),clim);syn;ylim([0 1e4]);xlim([-1.4 0.2])

%%%%%%%%%%
%%%%%%%%%%
% pu67bk2
%%%%%%%%%%
%%%%%%%%%%




%%%%%
% 10.27.10 - Summary figures
        load /cardinal4/r37g7/LearningCurves2.mat
        figure;hold on;
        plot(TClearnA1,1-e1.p05(1:end-1),'.')
        plot(TClearnA1,1-e1.p95(1:end-1),'.')
        plot(TClearnA1,1-e1.pmid(1:end-1),'.')
        for i=1:1000
            randsam=ceil(45*rand(1,45));
            guess(i)=mean(LClearn_apvA1(randsam)=='a');
        end
        plot([mean(TClearn_apvA1) mean(TClearn_apvA1)],[prctile(guess,5) prctile(guess,95)],'-','Color','r')
        plot(mean(TClearn_apvA1),mean(LClearn_apvA1=='a'),'.','Markersize',15,'Color','r')
        %
        plot(mean(TCbase_apv1)+460,mean(LCbase_apv1=='a'),'.','Markersize',15,'Color','r')
        guess=[];
        for i=1:1000
            randsam=ceil(length(LCbase_apv1)*rand(1,length(LCbase_apv1)));
            guess(i)=mean(LCbase_apv1(randsam)=='a');
        end
        plot([mean(TCbase_apv1)+460 mean(TCbase_apv1)+460],[prctile(guess,5) prctile(guess,95)],'-','Color','r')
        plot(mean(TCbase5)+560,mean(LCbase5=='a'),'.','Markersize',15,'Color','k')
        guess=[];
        for i=1:1000
            randsam=ceil(length(LCbase5)*rand(1,length(LCbase5)));
            guess(i)=mean(LCbase5(randsam)=='a');
        end
        plot([mean(TCbase5)+560 mean(TCbase5)+560],[prctile(guess,5) prctile(guess,95)],'-','Color','k')
        plot(mean(TCbase)+20,mean(LCbase=='a'),'.','Markersize',15,'Color','k')
        guess=[];
        for i=1:1000
            randsam=ceil(length(LCbase)*rand(1,length(LCbase5)));
            guess(i)=mean(LCbase(randsam)=='a');
        end
        plot([mean(TCbase)+20 mean(TCbase)+20],[prctile(guess,5) prctile(guess,95)],'-','Color','k')


% SyntaxFF - 11.01.10
%
% r37g7
%load /cardinal4/r37g7/Data428withFF.mat

% pu67bk2
datas=FFsyntax(3).data;
figure;hold on;
ravgwin=10;
for i=1:length(datas)
    if datas(i).acsf
        plot(runningmedian(timing3(datas(i).fvals),ravgwin),runningaverage(mean(datas(i).pitch(datas(1).window,:)),ravgwin),'*')
    else
        plot(runningmedian(timing3(datas(i).fvals),ravgwin),runningaverage(mean(datas(i).pitch(datas(1).window,:)),ravgwin),'r*')
    end
end
    


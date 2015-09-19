  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.18.2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cardinal4/crosscoINAs.mat
%%% 416.m
%%% figures in home directory

    %%% BENGALESE - new corr method (crossco)
        %%%% Bengalese Finch data - Lesions - get crossco data
            for i=1:9
                BFLresidsCTL(i).data=jc_residuals(BFLesion(i).UDpre20_1024);
            end
            for i=1:9
                BFLresidsINA(i).data=jc_residuals(BFLesion(i).UDpost20_1024);
            end
            BFLresidsCTL=BFLresidsCTL([1 3 4 9]);
            BFLresidsINA=BFLresidsINA([1 3 4 9]);
            BFLinitialpt=[800 400 320 240];
            BFLwidth=[350 450 300 200];
            for i=1:4
                clear nn
                for j=0:BFLwidth(i)
                    nn=corrcoef(BFLresidsINA(i).data(mean(BFLinitialpt(i+40)),:),BFLresidsINA(i).data(BFLinitialpt(i)+j,:));
                    crosscoBFLINA(i).data(j+1)=nn(2);
                end
            end
            for i=1:4
                clear nn
                for j=0:BFLwidth(i)
                    nn=corrcoef(BFLresidsCTL(i).data(BFLinitialpt(i),:),BFLresidsCTL(i).data(BFLinitialpt(i)+j,:));
                    crosscoBFLCTL(i).data(j+1)=nn(2);
                end
            end

        %%%% Bengalese Finch data - Inactivations that lower CV - get crossco data
        
        %%%%****%%% need to recalculate BFresid for inactivations unless I load the original data
            BFinitialpt=[450 450 300 700 200];
            BFwidth=[300 300 400 400 200];
            for i=1:5
                clear nn
                for j=0:BFwidth(i)
                    nn=corrcoef(BFresidINA(i).data(BFinitialpt(i),:),BFresidINA(i).data(BFinitialpt(i)+j,:));
                    crosscoBFINA(i).data(j+1)=nn(2);
                end
            end
            for i=1:5
                clear nn
                for j=0:BFwidth(i)
                    nn=corrcoef(BFresidCTL(i).data(BFinitialpt(i),:),BFresidCTL(i).data(BFinitialpt(i)+j,:));
                    crosscoBFCTL(i).data(j+1)=nn(2);
                end
            end
        %%%%% Calculate corr metric  - redundant
            %BFALLresidsCTL=[BFresidCTL BFLresidsCTL];

            BFALLinitialpt=[450 450 300 700 200 110 730 350 320 200];
            BFALLwidth=[300 300 400 400 200 200 350 400 250 200];
            BFALLinitialpt=BFALLinitialpt+50;

            for i=1:10
                clear nn
                for j=0:BFALLwidth(i)
                    nn=corrcoef(mean(BFALLresidsCTL(i).data(BFALLinitialpt(i):BFALLinitialpt(i)+39,:)),BFALLresidsCTL(i).data(BFALLinitialpt(i)+40+j,:));
                    crosscoBFCTL(i).data(j+1)=nn(2);
                end
            end
            for i=1:10
                clear nn
                for j=0:BFALLwidth(i)
                    nn=corrcoef(mean(BFALLresidsINA(i).data(BFALLinitialpt(i):BFALLinitialpt(i)+39,:)),BFALLresidsINA(i).data(BFALLinitialpt(i)+40+j,:));
                    crosscoBFINA(i).data(j+1)=nn(2);
                end
            end
        %%%% Plot BF
            figure;hold on;
            for i=1:10
                subplot(3,4,i,'XTickLabel',{'0','10','20','30'},'XTick',[0 80 160 240]);hold on;
                plot(crosscoBFCTL(i).data,'b')
                plot(crosscoBFINA(i).data,'r')
                xlim([0 320])
                ylim([0 1])
            end

   %%% ZEBRA - new corr method (crossco)         
        %%%% Zebra Finch Lesion data - get crossco data
            %%%% Get residuals %%
                for i=1:11
                    ZFresidsCTL(i).data=jc_residuals(ZFLesion(i).UDpre20_1024);
                end
                for i=1:11
                    ZFresidsINA(i).data=jc_residuals(ZFLesion(i).UDpost20_1024);
                end
            %%%%% Calculate corr metric %% 
                ZFinitialpt=[200 700 500 600 150 480 550 700 200 500 600];
                ZFwidth=[250 300 300 450 300 250 250 450 300 300 500];
                ZFinitialpt=ZFinitialpt+50;
                for i=1:11
                    clear nn
                    for j=0:ZFwidth(i)
                        nn=corrcoef(mean(ZFresidsCTL(i).data(ZFinitialpt(i):ZFinitialpt(i)+39,:)),ZFresidsCTL(i).data(ZFinitialpt(i)+40+j,:));
                        crosscoZFCTL(i).data(j+1)=nn(2);
                    end
                end
                for i=1:11
                    clear nn
                    for j=0:ZFwidth(i)
                        nn=corrcoef(mean(ZFresidsINA(i).data(ZFinitialpt(i):ZFinitialpt(i)+39,:)),ZFresidsINA(i).data(ZFinitialpt(i)+40+j,:));
                        crosscoZFINA(i).data(j+1)=nn(2);
                    end
                end
            %%%% Plot ZF
                figure;hold on;
                for i=1:11
                    subplot(3,4,i,'XTickLabel',{'0','10','20','30'},'XTick',[0 80 160 240]);hold on;
                    plot(crosscoZFCTL(i).data,'b')
                    plot(crosscoZFINA(i).data,'r')
                    ylim([0 1])
                    xlim([0 320])
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AUTOCORR PLOT - old method - middle 200 pts (25ms) %%%%%%%%%%
    % ZEBRA FINCHES
            Zinit=[200 720 520 700 200 500 580 800 250 520 750];
            Zfin=[400 920 720 900 400 700 780 1000 450 720 950];
            figure;subplot(221);hold on
            for j=1:11
                clear acC
                clear acI
                for i=1:20
                    acC(i,:)=xcorr(ZFresidsCTL(j).data(Zinit(j):Zfin(j),i));
                    acI(i,:)=xcorr(ZFresidsINA(j).data(Zinit(j):Zfin(j),i));
                end
                plot(mean(abs(acC))./max(mean(abs(acC))))
                XCwidthCTL(j)=(200-min(find(mean(abs(acC))./max(mean(abs(acC)))>0.5)))/8;
                plot(mean(abs(acI))./max(mean(abs(acI))),'r')
                XCwidthINA(j)=(200-min(find(mean(abs(acI))./max(mean(abs(acI)))>0.5)))/8;
            end
            figure;plot([XCwidthCTL;XCwidthINA],'-','Color','k')
            hold on;plot(1,XCwidthCTL,'*','Color','b')
            hold on;plot(2,XCwidthINA,'*','Color','r')
            %%%%%%
   % BENGALESE FINCHES  
            Binit=[500 550 350 850 145 140 850 400 350 200];
            Bfin=[700 750 550 1050 345 340 1050 600 550 400];

            subplot(222);hold on
            for j=1:10
                clear acC
                clear acI
                for i=1:20
                    acC(i,:)=xcorr(BFALLresidsCTL(j).data(Binit(j):Bfin(j),i));
                    acI(i,:)=xcorr(BFALLresidsINA(j).data(Binit(j):Bfin(j),i));
                end
                plot(mean(abs(acC))./max(mean(abs(acC))))
                BF_XCwidthCTL(j)=(200-min(find(mean(abs(acC))./max(mean(abs(acC)))>0.5)))/8;
                plot(mean(abs(acI))./max(mean(abs(acI))),'r')
                BF_XCwidthINA(j)=(200-min(find(mean(abs(acI))./max(mean(abs(acI)))>0.5)))/8;
            end

            subplot(223);hold on;
            plot([XCwidthCTL;XCwidthINA],'-','Color','k')
            plot(1,XCwidthCTL,'*','Color','b')
            plot(2,XCwidthINA,'*','Color','r')
            ylim([0 15])
            subplot(224);hold on;
            plot([BF_XCwidthCTL;BF_XCwidthINA],'-','Color','k')
            plot(1,BF_XCwidthCTL,'*','Color','b')
            plot(2,BF_XCwidthINA,'*','Color','r')
            ylim([0 15])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.18.2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at shape of learning - ChristmasFinal dataset
% program mat418.m
% saved in 418work.mat, figures in home directory

%%%% Final plot runs commands in 2 and 4
        Alldata1sig=[Alldatafirstten1sig Alldatasecondseven1sig];
        % 1. Baseline with learning overlayed
            figure;hold on;
            for i=1:17
                [pcurves]=mat418(Alldata1sig(i));
                [a,b]=hist(tfinal(i).data);
                a=(a./max(a))*200;
                subplot(4,5,i),'XTickLabel',{'0','20','40','60','80','100'},'XTick',[0 160 320 480 640 800]);hold on;
                xlim([1 1100]);
                plot(b,a+median(median(Alldata1sig(i).baselineAC)));
                plot(mean(Alldata1sig(i).baselineAC'),'k')
                plot(median(Alldata1sig(i).baselineAC'),'k')
                plot(mean(pcurves'),'r')
                plot(median(pcurves'),'r')
            end
        % 2. normalized - deltaPitch - plot this
            figure;hold on;
            for i=1:17
                [pcurves]=mat418(Alldata1sig(i));
                [a,b]=hist(tfinal(i).data);
                a=(a./max(a))*100;
                if strcmp(Alldata1sig(i).exp(1).direction,'down')
                    a=-1*a;
                end
                subplot(4,5,i);hold on%,'XTickLabel',{'0','20','40','60','80','100'},'XTick',[0 160 320 480 640 800]);hold on;
                xlim([1 1100]);
                plot(b,a);
                plot(mean(pcurves')-mean(Alldata1sig(i).baselineAC'),'r')
                plot(median(pcurves')-median(Alldata1sig(i).baselineAC'),'r')
            end
        % 3. How precise is variability? Autocorr function over middle 25ms
        
                xcwindowbegin=[300 270 340 270 230 300 500 450 500 240 230 400 550 600 950 800];
                % Get residuals
                for i=1:17
                    BFresids(i).data=jc_residuals(Alldata1sig(i).baselineAC);
                end
                figure;hold on
                    for j=1:17
                        clear acC
                        for i=1:50
                            acC(i,:)=xcorr(BFresids(j).data(xcwindowbegin(j):xcwindowbegin(j)+200,i));
                        end
                        plot(mean(abs(acC))./max(mean(abs(acC))))
                        BF_XCwidthCTL(j)=(200-min(find(mean(abs(acC))./max(mean(abs(acC)))>0.5)))/8;
                    end
        % 4. How precise is variability? crossco metric - plot this 
                % Get residuals
                for i=1:17
                    BFresids(i).data=jc_residuals(Alldata1sig(i).baselineAC);
                end
                for i=1:17
                    BFinitpt(i)=round(median(tfinal(i).data));
                end
                    BFALLwidth=[300 300 400 400 200 200 350 400 250 200];
                    for i=1:17
                        clear nn
                        for j=1:1700
                            nn=corrcoef(mean(BFresids(i).data(BFinitpt(i):BFinitpt(i)+39,:)),BFresids(i).data(j,:));
                            crosscoBFCTL(i).data(j+1)=nn(2);
                        end
                    end 
                hold on;
                %maxes used to normalize - get the real local extrema
                maxes=[165 320 -210 -300 330 -210 90 -115 115 175 230 85 -95 130 140 -100 135];
                for i=1:17
                    subplot(4,5,i);hold on;plot(crosscoBFCTL(i).data*maxes(i),'k');
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

% 42109
% ContingSim @70%
% Run #2 first 
for i=1:17
    ud=Alldata1sig(i).exp(1).direction;
    if strcmp(ud,'up')
        contingsim(i).data=ContingSim(BFinitpt(i),Alldata1sig(i).baselineAC,70);
    else
        contingsim(i).data=ContingSim(BFinitpt(i),Alldata1sig(i).baselineAC,30);
    end
    subplot(4,5,i);hold on;plot((contingsim(i).data-mean(Alldata1sig(i).baselineAC'))*(maxes(i)/maxesContingSim(i)),'k')
end
%%%%
%%% Asymmetry?
Bmiddle=[500 350 500 350 330 500 620 590 650 320 340 600 750 800 1100 950 580];
figure;
for i=1:17
    clear nn
    for j=1:1700
        nn=corrcoef(mean(BFresids(i).data(Bmiddle(i)-20:Bmiddle(i)+20,:)),BFresids(i).data(j,:));
        crosscoBFMiddle(i).data(j+1)=nn(2);
    end
    subplot(4,5,i);plot(crosscoBFMiddle(i).data);
end
figure;
for i=1:17
    clear nn
    for j=1:1700
        nn=corrcoef(mean(BFresids(i).data(Bmiddle(i)-20:Bmiddle(i)+20,:)),BFresids(i).data(j,:));
        crosscoBFMiddle(i).data(j+1)=nn(2);
    end
    subplot(4,5,i);plot(crosscoBFMiddle(i).data(Bmiddle(i)-200:Bmiddle(i)+200));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 4/23/09
%%%%%%%%%%%%%%%%%%%%%%%
% Predictability
% Predictability by INA vs. ALL

%%% Look at just the "nice" part of the note
for i=1:17
    s1(i)=Alldata1sig(i).startnote;
    s2(i)=Alldata1sig(i).endnote;
end
            figure;hold on;
            for i=1:17
                [pcurves]=mat418(Alldata1sig(i));
                [a,b]=hist(tfinal(i).data);
                a=(a./max(a))*100;
                if strcmp(Alldata1sig(i).exp(1).direction,'down')
                    a=-1*a;
                end
                subplot(4,5,i);hold on
                xlim([1 600]);
                plot(b-s1(i),a);
                absshiftedcurve(i).data=mean(pcurves(s1(i):s2(i),:)')-mean(Alldata1sig(i).baselineAC(s1(i):s2(i),:)');
                plot(mean(pcurves(s1(i):s2(i),:)')-mean(Alldata1sig(i).baselineAC(s1(i):s2(i),:)'),'r')
                plot(median(pcurves(s1(i):s2(i),:)')-median(Alldata1sig(i).baselineAC(s1(i):s2(i),:)'),'r')
            end
        % 4. How precise is variability? crossco metric - plot this 
                % Get residuals

                %maxes used to normalize - get the real local extrema
                maxes=[165 320 -210 -300 330 -210 90 -115 115 175 230 85 -95 130 140 -100 135];
                sign=[1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 1 1 -1 1];
                figure;hold on;
                for i=1:17
                    subplot(4,5,i);hold on
                    normshiftedcurve(i).data=(absshiftedcurve(i).data)./(sign(i)*max(abs((absshiftedcurve(i).data))));
                    plot((normshiftedcurve(i).data),'r')
                    plot(crosscoBFCTL(i).data(s1(i):s2(i)),'k');
                end
                

                for j=1:17
                distance=0;
                for k=1%0.8:0.01:1.2  % different heights of predictor function
                    long=length(normshiftedcurve(j).data);
                        distance=sum(abs(k*crosscoBFCTL(j).data(s1(j):s2(j))-normshiftedcurve(j).data));
                    distance=distance/long;
                end
                [dist(j),ind(j)]=min(distance);
                end
                
                
                    for i=1:17
                        clear nn
                        for j=1:1700
                            nn=corrcoef(mean(BFresids(i).data(BFinitpt(i)-s1p(i):BFinitpt(i)-s1p(i)+39,:)),BFresids(i).data(j,:));
                            crosscoBFCTL(i).data(j+1)=nn(2);
                        end
                    end 

%%%%% 4/24/09
figure;
for i=1:12
    subplot(3,4,i);hold on;plot(BFresids(Longnotes(i)).data)
    plot(tfinal(Longnotes(i)).data,0.2,'*','Color','k')
    plot(std(BFresids(Longnotes(i)).data'),'k','LineWidth',4)
end
for j=1:1400
        nn=corrcoef(mean(g5o78.baselineAC(880-20:880+20,:)),g5o78.baselineAC(j,:));
        crossco(j+1)=nn(2);
end

%%%%% 4/26/09
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(toffs1)
    for j=1:1400
            nn=corrcoef(ResidB61inas(round(toffs1(k)),:),ResidB61inas(j,:));
            crosscoB(k,j+1)=nn(2);
    end
end
crosscorrB=mean(crosscoB);


for k=1:length(toffsbk63)
    for j=1:1400
            nn=corrcoef(ResidB63ina(round(toffsbk63(k)),:),ResidB63ina(j,:));
            crossco63B(k,j+1)=nn(2);
    end
end
crosscorr63B=mean(crossco63B);

for k=1:length(toffsbk63)
    for j=1:1400
            nn=corrcoef(ResidB63(round(toffsbk63(k)),:),ResidB63(j,:));
            crossco63A(k,j+1)=nn(2);
    end
end
crosscorr63A=mean(crossco63A);


% match metric
% scaling
sides=[364 980]; % liberal

for i=1:50
    fac=100+2*i;
    matchscoreB(i)=sum(abs(crosscorrB(sides(1):sides(2))-learningBK61(sides(1):sides(2))./fac));
end
for i=1:50
    fac=100+2*i;
    matchscoreA(i)=sum(abs(crosscorrA(sides(1):sides(2))-learningBK61(sides(1):sides(2))./fac));
end

sidesBK63=[380 850]; % liberal
for i=1:50
    fac=150+2*i;
    matchscore63A(i)=sum(abs(crosscorr63A(sidesBK63(1):sidesBK63(2))-learningBK63(sidesBK63(1):sidesBK63(2))./fac));
end
for i=1:50
    fac=150+2*i;
    matchscore63B(i)=sum(abs(crosscorr63B(sidesBK63(1):sidesBK63(2))-learningBK63(sidesBK63(1):sidesBK63(2))./fac));
end
figure;plot(matchscoreA);hold on;plot(matchscoreB,'r')
[mA,indA]=min(matchscoreA);
[mB,indB]=min(matchscoreB);
[m63A,ind63A]=min(matchscore63A);
[m63B,ind63B]=min(matchscore63B);
facA=indA*2+100;
facB=indB*2+100;
fac63A=ind63A*2+150;
fac63B=ind63B*2+150;
figure;hold on;
plot(crosscorr63A*fac63A*-1);plot(crosscorr63B*fac63B*-1,'r');plot(learningBK63,'k')


% pu34
for k=1:length(toffspu34)
    for j=1:1400
            nn=corrcoef(ResidPU34ina(round(toffspu34(k)),:),ResidPU34ina(j,:));
            crosscoPU34ina(k,j+1)=nn(2);
    end
end
crosscorrPU34ina=mean(crosscoPU34ina);
for k=1:length(toffspu34)
    for j=1:1400
            nn=corrcoef(ResidPU34acsf(round(toffspu34(k)),:),ResidPU34acsf(j,:));
            crosscoPU34ac(k,j+1)=nn(2);
    end
end
crosscorrPU34ac=mean(crosscoPU34ac);

for k=1:length(PU34mockB63toffs)
    for j=1:1400
            nn=corrcoef(ResidBaselineB63(round(PU34mockB63toffs(k)),:),ResidBaselineB63(j,:));
            crosscoB63mPU34ac(k,j+1)=nn(2);
    end
end
crosscorrB63mPU34ac=mean(crosscoB63mPU34ac);


%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%4/27/09

for i=1:17
    [pcurves]=mat418(Alldata1sig(i));
    learned(i).data=mean(pcurves')-mean(Alldata1sig(i).baselineAC');
end
for i=1:26
    %%%% Get self-predictions
    clear crossco;
    for k=1:100
        Targetrandomdraw=round(rand*(length(Predict(i).Targeting)-1)+1);
        for j=1:1400
            nn=corrcoef(Predict(i).ResidAC(round(Predict(i).Targeting(Targetrandomdraw)),:),Predict(i).ResidAC(j,:));
            crossco(k,j+1)=nn(2);
        end
    end
    Predict(i).crosscorr=mean(crossco);
    subplot(5,6,i);plot(Predict(i).crosscorrAC);hold on;plot(Predict(i).Learned./abs(maxes(Longnotes(i))))
end
    % optimize normalization factor
    % ****** Do this again with better (more conservative) s1 and s2 ---
    % ****plus add in bk63 and bk61
    for j=1:14
         for i=1:50
            fac=maxes(Longnotes(j))-30+2*i;
            matchscore(i)=sum(abs(Predict(j).crosscorr(s1(Longnotes(j)):s2(Longnotes(j)))-Predict(j).Learned(s1(Longnotes(j)):s2(Longnotes(j)))./fac));
         end
         [m,ind]=min(matchscore);
         facF(j)=ind*2-30+maxes(Longnotes(j));
         subplot(4,4,j);plot(Predict(j).crosscorr);hold on;plot(Predict(j).Learned./abs(facF(j)))
    end
    %%%%%%%%%
    
    
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%% 4.28.09     %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%% See cardinal/428.mat

%bkw's
for i=9
    %%%% Get self-predictions
    clear crossco;
    for k=1:100
        Targetrandomdraw=round(rand*(length(Predict(i).Targeting)-1)+1);
        for j=1:1400
            nn=corrcoef(Predict(i).ResidAC(round(Predict(i).Targeting(Targetrandomdraw)),:),Predict(i).ResidAC(j,:));
            crossco(k,j+1)=nn(2);
        end
    end
    Predict(i).crosscorrAC=mean(crossco);
end
%%%
figure;
for i=1:14
subplot(4,4,i);plot(Predict(i).crosscorrAC(Predict(i).onset:Predict(i).offset),'r');hold on;plot(Predict(i).Learned(Predict(i).onset:Predict(i).offset)./abs(mmax(i)),'k')
end  
%%%%%%%%%%%
%%% Optimize prediction
%%%%%%%%%%%
figure;
    for j=1:14
         for i=1:200
            fac=abs(mmax(j))-60+2*i;
            matchscore(i)=sum(abs(Predict(j).crosscorrAC(Predict(j).onset:Predict(j).offset)-Predict(j).Learned(Predict(j).onset:Predict(j).offset)./fac));
         end
         [best(j),ind]=min(matchscore);
         facF(j)=ind*2-60+abs(mmax(j));
         Predict(j).LearnedNorm=Predict(j).Learned./abs(facF(j));
         subplot(4,4,j);plot(Predict(j).crosscorrAC(Predict(j).onset:Predict(j).offset),'r');hold on;plot(Predict(j).LearnedNorm(Predict(j).onset:Predict(j).offset))
    end
    %%%%%%%%%
%%%% R2 metric - how good is prediction?
for j=1:14
    Observed=Predict(j).LearnedNorm(Predict(j).onset:Predict(j).offset);
    Expected=Predict(j).crosscorrAC(Predict(j).onset:Predict(j).offset);
    r2(j)=1-sum((Observed-Expected).^2)./sum((Observed-mean(Observed)).^2);
end
%%%%%%%%%
%%%%%%%%%
% compare with or92 (bestest)

% Get best for the right region of the original bird

%test
for j=[1:16 19]
    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
    for i=[1:16 19]
        NewTargs=Predict(i).onset+NormTargs;
        clear crossco;
        Prediction(i).data=jccrossco(Predict(i).ResidAC,NewTargs); % round?
        % best prediction
        default=1./max(Prediction(i).data(Predict(i).onset:Predict(i).offset));
        for x=1:200
            fac(x)=default-0.8+0.005*x;
            matchscore(x,:)=sum(abs(Prediction(i).data(Predict(i).onset:Predict(i).onset+350)*fac(x)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
        end
        [BestSim3(i,j),ind]=min(matchscore);
    end
end
%%%%% 
%%%%% Fitting to #12 (targ @ end) indicates:
%%%%% Fast timescale: (12,13,14) - or92or82,bk61w42,bk63w63
%%%%% Intermediate timescales: (4/5/6, 8/9) - bk50w18,pk37bk19
%%%%% Slow timescale: (1/2,3,7,10,11) - bk15bk14,bk20bk45,pk32,pk39,pu34
% INHERITED? TUTOR?
% To do - check across birds (i.e. not just 12).
% Clean up data (428.mat)
% Check LMAN
% Screen ZF song.

%%%%%%
%%% 4.29.09 %%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%% April 29, 2009
%%%%%%%%%%%%%%%%%%%

% Comparison between syllables

% 1. Check FWHM - agree with Fast/Slow evaluation?
figure;
for j=[1:16 17]
    startpt=Predict(j).onset;
    endpt=Predict(j).onset+300;
    clear ac
    for i=1:size(Predict(j).ResidAC,2)
    data(j).ac(i,:)=xcorr(Predict(j).ResidAC(startpt:endpt,i));
    end
%     gplot=median(abs(ac))./max(median(abs(ac)));
%     hold on;plot(gplot)
%     ab(j)=300-min(find(gplot>0.5));
end
figure;plot(ab)
% 2. Fit to #1- Targ at beginning
bird=1;
for j=bird
    for i=1:200
        fac=abs(mmax(j))-60+2*i;
        matchscore(i)=sum(abs(Predict(j).crosscorrAC(Predict(j).onset:Predict(j).onset+350)-Predict(j).Learned(Predict(j).onset:Predict(j).onset+350)./fac));
    end
    [best(bird),ind]=min(matchscore);
    facF(j)=ind*2-60+abs(mmax(j));
    Best1=Predict(j).Learned./abs(facF(j));
end

NormTargs=Predict(bird).Targeting-Predict(bird).onset; % dist b/f offset
%longness=Predict(bird).offset-Predict(bird).onset+1; %%% unnecessary - using 350
figure;hold on;plot(Best1(Predict(bird).onset:Predict(bird).onset+350),'k');plot(Predict(bird).crosscorrAC(Predict(bird).onset:Predict(bird).onset+350),'r')
for i=2:14
    NewTargs=Predict(i).onset+NormTargs;
    % init prediction
    clear crossco;
    for k=1:100
        Targetrandomdraw=round(rand*(length(NewTargs)-1)+1);
        for j=1:1400
            nn=corrcoef(Predict(i).ResidAC(round(NewTargs(Targetrandomdraw)),:),Predict(i).ResidAC(j,:));
            crossco(k,j+1)=nn(2);
        end
    end
    Prediction(i).data=mean(crossco);
    % best prediction
    for ii=1:200
        fac=abs(mmax(bird))-60+2*ii;
        matchscore(ii)=sum(abs(Prediction(i).data(Predict(i).onset:Predict(i).onset+350)-Predict(bird).Learned(Predict(bird).onset:Predict(bird).onset+350)./fac));
    end
    [best(i),ind]=min(matchscore);
    facF1=ind*2-60+abs(mmax(bird));
    Prediction(i).bestdata=Prediction(i).data*(facF1/facF(bird));
    plot(Prediction(i).bestdata(Predict(i).onset:Predict(i).onset+350),'b')
end

%%%%%%%%%%%%%%%%%%%
%%%% April 30, 2009
%%%%%%%%%%%%%%%%%%%

% LMAN only vs. ALL
INAs=[3 7 10 11 13 14:22];
figure;
for i=1:15
    subplot(4,4,i);hold on;
    plot(Predict(INAs(i)).LearnedNorm,'k')
    plot(Predict(INAs(i)).crosscoAC/max(Predict(INAs(i)).crosscoAC(Predict(INAs(i)).onset:Predict(INAs(i)).offset)),'b')
    plot(Predict(INAs(i)).crosscoINA/max(Predict(INAs(i)).crosscoINA(Predict(INAs(i)).onset:Predict(INAs(i)).offset)),'r')
    plot(Predict(INAs(i)).crosscorrAC-(Predict(INAs(i)).crosscorrINA-Predict(INAs(i)).crosscorrAC),'g')
    xlim([Predict(INAs(i)).onset Predict(INAs(i)).offset]);ylim([0 1.05])
end
%%%%%%%%%%%%%%%%%%
%%%% May 1, 2009
%%%%%%%%%%%%%%%%%%
numexps=3;
% LearnedNorm is normalized to one
% 1. Calculate crossco predictions
for i=[3 7 10 11 13 14]
    Predict(i).crosscoAC=jccrossco(Predict(i).ResidAC,Predict(i).Targeting);
    if ~isempty(Predict(i).ResidINA)
        Predict(i).crosscoINA=jccrossco(Predict(i).ResidINA,Predict(i).Targeting);
    end
end
% 2. Optimize crossco predictions
for i=1:22
    [Predict(i).bestcrosscorrAC,Predict(i).mindistAC]=optimalcc(Predict(i),1);
    if ~isempty(Predict(i).crosscoINA)
        [Predict(i).bestcrosscorrINA,Predict(i).mindistINA]=optimalcc(Predict(i),0);
    end
end
% Note that in 5/6 cases (3,7,11,13,14 but not 10), mindistAC < mindistINA
% 3. Plot predictions normalized to peak (not necessarily best ones)
figure;hold on;
a=0;
for i=[3 7 10 11 13:22]
    a=a+1;
    subplot(4,4,a);hold on;
    plot(Predict(i).LearnedNorm,'k')
    xlim([Predict(i).onset Predict(i).offset]); ylim([0 1.05])
    defaultAC=1./max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
    plot(Predict(i).crosscoAC*defaultAC)
    if ~isempty(Predict(i).crosscoINA)
        defaultINA=1./max(Predict(i).crosscoINA(Predict(i).onset:Predict(i).offset));
        plot(Predict(i).crosscoINA*defaultINA,'r')
        plot(2*Predict(i).crosscoAC*defaultAC-Predict(i).crosscoINA*defaultINA,'g')
    end
end


% May 4th - prep for lab meeting?
inas=[3 7 10 11 13:22];
for i=1:length(inas)
    AcFit(i)=Predict(inas(i)).mindistAC/(Predict(inas(i)).offset-Predict(inas(i)).onset+1);
    InaFit(i)=Predict(inas(i)).mindistINA/(Predict(inas(i)).offset-Predict(inas(i)).onset+1);
end
figure;plot([AcFit;InaFit],'-','Color','k')
hold on;plot(1,AcFit,'*','Color','b')
hold on;plot(2,InaFit,'*','Color','r')
hold on;plot([mean(AcFit),mean(InaFit)],'-','Color','k','Linewidth',4)
hold on;plot(2,mean(InaFit),'+','Color','k')
hold on;plot(1,mean(AcFit),'+','Color','k')
% Autocorrelation metric generally agrees.
inas=[7 13 14 15 16 19 24:28];
for ii=1:11
    k=inas(ii);
    clear ac
    clear ina
    for i=1:size(Predict(k).ResidAC,2)
        ac(i,:)=((xcorr(Predict(k).ResidAC(Predict(k).onset+100:Predict(k).onset+300,i))));
    end
    for i=1:size(Predict(k).ResidINA,2)
        ina(i,:)=((xcorr(Predict(k).ResidINA(Predict(k).onset+100:Predict(k).offset-100,i))));
    end
    ccac(k).data=mean(ac)./max(mean(ac));
    ccina(k).data=mean(ina)./max(mean(ina));
    dat(ii)=sum(ccac(k).data-ccina(k).data)/length(ccac(k).data);
end
%%%% ****** "abs(dat)" is highly correlated (0.6) with CV, especially for the
%%%% later ones that I've cleaned up the means on (0.8). ******* The
%%%% strength of this is the same whether you leave in or exclude the
%%%% negative valued ones, suggesting that it isn't just tighter contours
%%%% lowering CV.

%%%% Prediction from itself
    figure;hold on;
    a=0;
    for i=[3 7 10 11 13:22]
        a=a+1;
        subplot(4,4,a);hold on;
        plot(Predict(i).LearnedNorm,'k')
        xlim([Predict(i).onset Predict(i).offset]); ylim([0 1.05])
        defaultAC=1./max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
        plot(Predict(i).crosscoAC*defaultAC)
        if ~isempty(Predict(i).crosscoINA)
            defaultINA=1./max(Predict(i).crosscoINA(Predict(i).onset:Predict(i).offset));
            plot(Predict(i).crosscoINA*defaultINA,'r')
            plot(2*Predict(i).crosscoAC*defaultAC-Predict(i).crosscoINA*defaultINA,'g')
        end
    end

%%%% Compare to predictions from other ones
    bird=1; %%%% DO THIS FOR EACH
    NormTargs=Predict(bird).Targeting-Predict(bird).onset; % dist b/f offset
    for i=[3 7 10:14]
        NewTargs=Predict(i).onset+NormTargs;
        % init prediction
        clear crossco;
        for k=1:100
            Targetrandomdraw=round(rand*(length(NewTargs)-1)+1);
            for j=1:1400
                nn=corrcoef(Predict(i).ResidAC(round(NewTargs(Targetrandomdraw)),:),Predict(i).ResidAC(j,:));
                crossco(k,j+1)=nn(2);
            end
        end
        Prediction(i).bestdata=mean(crossco);
        % best predictbion
        plot(Prediction(i).bestdata(Predict(i).onset:Predict(i).onset+350),'b')
    end
    % plot
    figure;hold on;
    for i=[3 7 10:14]
        plot(Prediction(i).bestdata(Predict(i).onset:Predict(i).onset+350)/max(Prediction(i).bestdata(Predict(i).onset:Predict(i).onset+350)),'b')
    end
    hold on;plot(Predict(1).crosscoAC(Predict(1).onset:Predict(1).onset+350)/max(Predict(1).crosscoAC(Predict(1).onset:Predict(1).onset+350)),'r')
    hold on;plot(Predict(1).LearnedNorm(Predict(1).onset:Predict(1).onset+350),'k')
%%% It's worth considering only looking later in the note (((i.e. not right
%%% at the onset, because prediction tends to be crappier near the onset)))
inas=[7 11 13:23];
for ii=1:15
    ons=Predict(inas(ii)).onset;
    offs=Predict(inas(ii)).offset;
    CV(ii)=1-mean(std(Predict(inas(ii)).ResidINA(ons:offs,:)')./std(Predict(inas(ii)).ResidAC(ons:offs,:)'));
end
%%% CV vs abs(dat) plot
figure;subplot(221);plot(CV,(dat),'*')
hold on;plot([-1 1],[0 0],'-','Color','k')
hold on;plot(CV(find(dat<0)),abs(dat(find(dat<0))),'*','Color','r')
xlim([-0.1 0.6]);ylim([-0.05 0.06])
subplot(223);plot(CV,abs(dat),'*')
hold on;plot([-1 1],[0 0],'-','Color','k')
xlim([-0.1 0.6]);ylim([-0.05 0.06])
subplot(222);plot(CV,'*')
hold on;plot(find(dat<0),CV(find(dat<0)),'*','Color','r')

%%%%% Zebra Finch data set %%%%%%%%%
for i=1:11
    ZFLesion(i).ResidUDpre20=jc_residuals(ZFLesion(i).UDpre20_1024);
    ZFLesion(i).ResidUDpost20=jc_residuals(ZFLesion(i).UDpost20_1024);
end
for k=1:11
    clear pre
    clear post
    for i=1:size(ZFLesion(k).ResidUDpre20,2)
        pre(i,:)=((xcorr(ZFLesion(k).ResidUDpre20(ind1ZF(k)-130:ind1ZF(k)+130,i))));
    end
    for i=1:size(ZFLesion(k).ResidUDpost20,2)
        post(i,:)=((xcorr(ZFLesion(k).ResidUDpost20(ind1ZF(k)-130:ind1ZF(k)+130,i))));
    end
    ccpre(k).data=mean(pre)./max(mean(pre));
    ccpost(k).data=mean(post)./max(mean(post));
end
for k=1:11
    datLes(k)=sum(ccpre(k).data-ccpost(k).data)/length(ccpre(k).data);
end
% r^2=0.95!!!! between datLes and datDir - and all but one have broad DIR
% song and broad LES song
for k=1:6
    datDir(k)=sum(ccpre(k).data-ccdir(k).data)/length(ccpre(k).data);
end
for ii=1:11
    ons=ind1ZF-200;
    offs=ind1ZF+200;
    CV(ii)=1-mean(std(ZFLesion(ii).ResidUDpost20(ons:offs,:)')./std(ZFLesion(ii).ResidUDpre20(ons:offs,:)'));
end
% probable relationship between CV and var reduction for dir song (quite poor for lesion song)

%%%%% 5.11.09
% What...is LMAN doing?
t=-5:0.02:5;
count=0;
for i=[7 11 13:22];
    % Normalize crosscoAC and crosscoINA to peak
        ncrosscoAC=Predict(i).crosscoAC./max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
        ncrosscoINA=Predict(i).crosscoINA./max(Predict(i).crosscoINA(Predict(i).onset:Predict(i).offset));
    % try a bunch of different weightings of the difference
        for j=1:length(t)
            tt=t(j);
            ncrosscoMAN=ncrosscoAC+tt*(ncrosscoAC-ncrosscoINA);
            default=1;
            for k=1:200
                fac(k)=default-0.8+0.005*k;
                matchscore(k)=sum(abs(ncrosscoMAN(Predict(i).onset:Predict(i).offset-100)*fac(k)-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset-100)));
            end
            [mindist(j),ind]=min(matchscore);
        end
        [a,b]=min(mindist);
        count=count+1;
        bestfactor(count)=t(b);
    % for each weighting, optimize height
end

% Predictive ability from mutual?
    for j=[1:16 19]
    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
        for i=[1:16 19]
            NewTargs=Predict(i).onset+NormTargs;
            clear crossco;
            Prediction(i).data=jccrossco(Predict(i).ResidAC,NewTargs); % round?
            % best prediction
            default=1./max(Prediction(i).data(Predict(i).onset:Predict(i).offset));
            for x=1:200
                fac(x)=default-0.8+0.005*x;
                matchscore(x,:)=sum(abs(Prediction(i).data(Predict(i).onset:Predict(i).onset+350)*fac(x)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
            end
            [BestSim3(i,j),ind]=min(matchscore);
        end
    end
    g=[1:17];
    for i=1:17
        BestScoreSelf(i)=BestSims2(i,i);
        BestScoreOther(i)=mean(BestSims2(find(g~=i),i));
    end
    [h,p]=ttest(BestScoreSelf,BestScoreOther,[],'left')
    figure;plot(1./BestScoreOther,1./BestScoreSelf,'*')
    c=0:0.001:1;
    d=0:0.001:1;
    hold on;plot(c,d,'Color','k')
    xlim([0 0.16])
    ylim([0 0.16])
    hold on;plot(mean(1./BestScoreOther),mean(1./BestScoreSelf),'+')
 % Predictive ability w/o targeting information?
    for j=1:22
            NormTargs=median(Predict(j).Targeting); % dist b/f offset
            clear crossco;
            Prediction(j).data=jccrossco2(Predict(j).ResidAC,NormTargs); % round?
            % best prediction
            default=1./max(Prediction(j).data(Predict(j).onset:Predict(j).offset));
            for x=1:200
                fac(x)=default-0.8+0.005*x;
                matchscore(x,:)=sum(abs(Prediction(j).data(Predict(j).onset:Predict(j).onset+350)*fac(x)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
            end
            [BestSim4(j),ind]=min(matchscore);
    end
    figure;plot(1./BestScoreNoTarg,1./BestScoreSelf,'*')
    c=[0:0.1:50];
    d=c;
    hold on;plot(c,d,'Color','k')
    xlim([0 50])
    ylim([0 50])
    hold on;plot(mean(1./BestScoreNoTarg),mean(1./BestScoreSelf),'+')
 % Predictive ability of DC shift?
    for i=[1:16 19]
        mindistDC(i)=sum(abs(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350)-mean(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350))));
    end
    figure;plot(1./BestScoreDC,1./BestScoreSelf,'*')
    hold on;plot(c,d,'Color','k')
    hold on;plot(mean(1./BestScoreDC),mean(1./BestScoreSelf),'+')
 % Predictive ability from LMAN lesioned/inactivated
figure;plot(1./BestScoreLESION,1./BestScoreINTACT,'*')
hold on;plot(c,d,'Color','k')
ylim([0 20])
xlim([0 20])
hold on;plot(mean(1./BestScoreLESION([1:6 9])),mean(1./BestScoreINTACT([1:6 9])),'+')
[h,p]=ttest(BestScoreINTACT,BestScoreLESION,[],'left')
%%%% Plot all
    figure;hold on;
    a=0;
    for i=1:26
        a=a+1;
        subplot(5,6,a);hold on;
        plot(Predict(i).LearnedNorm,'k')
        xlim([Predict(i).onset Predict(i).offset]); ylim([0 1.05])
        defaultAC=1./max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
        plot(Predict(i).crosscoAC*defaultAC)
        if ~isempty(Predict(i).crosscoINA)
            defaultINA=1./max(Predict(i).crosscoINA(Predict(i).onset:Predict(i).offset));
            plot(Predict(i).crosscoINA*defaultINA,'r')
        end
    end
    
% median targ and max shift relationship
    xx=[0:1:100];
    for i=[1:22]
        medtarg(i)=median(Predict(i).Targeting)-Predict(i).onset;
        [a,b]=max(abs(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset)));
        maxshift(i)=b;
    end
    for i=1:22
        stdtargs(i)=std(Predict(i).Targeting);
    end
    figure;errorbar(medtarg/8,maxshift/8,stdtargs/8,'*')
    hold on;plot(xx,xx,'-','Color','k')
    xlim([0 70]); ylim([0 70])
    % r=0.78

% as proportion of note
    xx=[0:1:100];
    for i=1:22
        ll=Predict(i).offset-Predict(i).onset;
        medtarg(i)=(median(Predict(i).Targeting)-Predict(i).onset)/ll;
        [a,b]=max(abs(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset)));
        maxshift(i)=b/ll;
    end
    % only slightly worse - r=0.7
% Center of mass relationship
    for i=1:22
        cmasstarg(i)=median(Predict(i).Targeting)-Predict(i).onset;
        ll=Predict(i).offset-Predict(i).onset;
        rundist=[];
        for j=1:ll
            rundist(j)=Predict(i).LearnedNorm(Predict(i).onset+j)*j;
        end
        cmassshift(i)=sum(rundist)/sum(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    end
    % This is worse, but still barely significant - r=0.5
%%% WIDTH METRIC
% for i=1:22
%     % 20ms
%     prednorm=Prediction(i).data/max(Prediction(i).data(Predict(i).onset:Predict(i).offset));
%     [a,maxim]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
%     if maxim-160>0
%         decayLeftActual(i)=1-Predict(i).LearnedNorm(Predict(i).onset+maxim-160);
%     end
%     if maxim+160<Predict(i).offset-Predict(i).onset
%         decayRightActual(i)=1-Predict(i).LearnedNorm(Predict(i).onset+maxim+160);
%     end
%     [a,maxPredict]=max(prednorm(Predict(i).onset:Predict(i).offset));
%     if maxPredict-160>0
%         decayLeftPredict(i)=1-prednorm(Predict(i).onset+maxPredict-160);
%     end
%     if maxPredict+160<Predict(i).offset-Predict(i).onset
%         decayRightPredict(i)=1-prednorm(Predict(i).onset+maxPredict+160);
%     end
% end
%%%% QUANTIFY temporal restrictedness of shift
decays=cc1c(Predict);
figure;hold on;
for i=1:length(decays)
    decayed(i)=decays(i).data(end);
end
%%%% QUANTIFY temporal restrictedness of variability - predict based on
%%%% targeting
% Prediction is the prediction w/o targeting distn
decayCorr=cc2c(Predict,Prediction);
for i=1:22
decayedCorr(i)=decayCorr(i).data(end);
end
figure;plot(1,decayedCorr([1:14 18:22]),'*','Color','k')
hold on;plot(2,decayed([1:14 18:22]),'*','Color','k')
hold on;plot([decayedCorr([1:14 18:22]);decayed([1:14 18:22])],'Color','k')
hold on;plot(1,mean(decayedCorr([1:14 18:22])),'+')
hold on;plot(2,mean(decayed([1:14 18:22])),'+')
hold on;plot([mean(decayedCorr([1:14 18:22]));mean(decayed([1:14 18:22]))])

figure;plot(1,1./BestScoreDC,'*','Color','k')
hold on;plot(1,mean(1./BestScoreDC),'+','Color','k')
hold on;plot(2,1./BestScoreNoTarg,'*','Color','k')
hold on;plot(2,mean(1./BestScoreNoTarg),'+','Color','k')
hold on;plot(3,1./BestScoreSelf,'*','Color','k')
hold on;plot(3,mean(1./BestScoreSelf),'+','Color','k')
hold on;plot(4,1./BestScoreOther,'*','Color','k')
hold on;plot(4,mean(1./BestScoreOther),'+','Color','k')


%%%% GENERALLY AGREES WITH CONTING SIM
a=0
for i=1:22
    a=a+1;subplot(5,5,a);hold on;
    CSs(i,j).data=(ContingSim2(Predict(i).Targeting,Predict(i).ResidAC,70));
figure;hold on;
for i=1:26
    subplot(5,6,i)
    aa=CSs(i,j).data; 
    plot(aa/max(aa(Predict(i).onset:Predict(i).offset)),'r')
    hold on;plot(Predict(i).LearnedNorm,'k')
    hold on;plot(Predict(i).bestcrosscorrAC./max(Predict(i).bestcrosscorrAC(Predict(i).onset:Predict(i).offset)))
    xlim([Predict(i).onset Predict(i).offset]);
end

%%%% Mutual predictions
for i=1:17
    for j=1:17
    aa=corrcoef(BestSims2(i,:),BestSims2(j,:));
    mutual(i,j)=aa(2);
    end
end
%%%% 2,16,14 vs. 10,12,15,17
k=0;
for i=1:17
    k=k+1;
    midd=200;
    atest(:,k)=jccrossco2(Predict(i).ResidAC(Predict(i).onset:Predict(i).onset+400,:),midd);
end

clear mutual2
ii=0;
for i=find(find(ab>55)<65)
    ii=ii+1;
    jj=0;
    for j=find(find(ab>55)<65)
        jj=jj+1;
        aa=corrcoef(BestSims2(i,:),BestSims2(j,:));
        mutual2(ii,jj)=aa(2);
    end
end
mean(mean(mutual2))

%%%% 5.20.09
fronts=[1:3 5:7 10:11 13:15 18:19 21:26];
ends=[4 8 9 12 16 17 20];
for i=[1:3 5:7 10:11 13:15 18:19 20 22]
    [maxpred,predtop]=max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
    [maxact,acttop]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    dropPredict(i)=Predict(i).crosscoAC(predtop+Predict(i).onset+200)/maxpred;
    dropActual(i)=Predict(i).LearnedNorm(acttop+Predict(i).onset+200)/maxact;
end
for i=[4 8 9 12 16 17 21]
    [maxpred,predtop]=max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
    [maxact,acttop]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    dropPredict(i)=Predict(i).crosscoAC(predtop+Predict(i).onset-200)/maxpred;
    dropActual(i)=Predict(i).LearnedNorm(acttop+Predict(i).onset-200)/maxact;
end
% Short notes
for i=[23:26]
    [maxpred,predtop]=max(Predict(i).crosscoAC(Predict(i).onset:Predict(i).offset));
    [maxact,acttop]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    dropPredict(i)=Predict(i).crosscoAC(predtop+Predict(i).onset+120)/maxpred;
    dropActual(i)=Predict(i).LearnedNorm(acttop+Predict(i).onset+120)/maxact;
end
figure;subplot(212);plot(dropPredict,dropActual,'*')
hold on;plot(dropPredict(23:26),dropActual(23:26),'*','Color','r')
xlim([0.45 0.95])
ylim([0.45 0.95])
subplot(211);hold on;plot(dropPredict,dropActual,'*')
hold on;plot(dropPredict(23:26),dropActual(23:26),'*','Color','r')

%%% 5.21.09
for i=1:26
    Predict(i).crosscoACnotarg=jccrossco2(Predict(i).ResidAC,mean(Predict(i).Targeting));
end
for i=[1:3 5:7 10:11 13:15 18:19 20 22]
    [maxpred,predtop]=max(Predict(i).crosscoACnotarg(Predict(i).onset:Predict(i).offset));
    dropPredictnt(i)=Predict(i).crosscoACnotarg(predtop+Predict(i).onset+200)/maxpred;
end
for i=[4 8 9 12 16 17 21]
    [maxpred,predtop]=max(Predict(i).crosscoACnotarg(Predict(i).onset:Predict(i).offset));
    dropPredictnt(i)=Predict(i).crosscoACnotarg(predtop+Predict(i).onset-200)/maxpred;
end
% Short notes
for i=[23:26]
    [maxpred,predtop]=max(Predict(i).crosscoACnotarg(Predict(i).onset:Predict(i).offset));
    dropPredictnt(i)=Predict(i).crosscoACnotarg(predtop+Predict(i).onset+120)/maxpred;
end

% Predictive ability from mutual?
% BestSimX(row 2 tells how well syllable 2 predicts all other notes
%         (column 2 tells how well syllable 2 is predicted by all notes)
    for j=1:22
    %NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
        for i=1:22
            NewTargs=Predict(i).onset+NormTargs;
            %clear crossco;
            %Prediction(i,j).data=jccrossco(Predict(i).ResidAC,NewTargs(find(NewTargs>64))); % round?
            % best prediction
            default=1./max(Prediction(i,j).data(Predict(i).onset:Predict(i).offset));
            for x=1:200
                fac(x)=default-0.8+0.005*x;
                matchscore(x,:)=sum(abs(Prediction(i,j).data(Predict(i).onset:Predict(i).onset+350)*fac(x)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
            end
            [BestSimX(i,j),ind]=min(matchscore);
        end
    end
    for j=1:22
    %NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
        for i=1:22
            %NewTargs=Predict(i).onset+NormTargs;
            %clear crossco;
            %Prediction(i,j).data=jccrossco(Predict(i).ResidAC,NewTargs(find(NewTargs>64))); % round?
            % best prediction
            default(i,j)=1./max(Prediction(i,j).data(Predict(i).onset:Predict(i).offset));
                matchscore=sum(abs(Prediction(i,j).data(Predict(i).onset:Predict(i).onset+350)*default(i,j)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));

            [BestSimX(i,j)]=(matchscore);
        end
    end
    
    g=[1:22];
    for i=1:22
        BestScoreSelf(i)=mean(BestSimX(i,i));
        BestScoreOther(i)=mean(BestSimX(find(g~=i),i));
     end
%     [h,p]=ttest(BestScoreSelf,BestScoreOther,[],'left')
%     figure;plot(1./BestScoreOther,1./BestScoreSelf,'*')
%     c=0:0.001:1;
%     d=0:0.001:1;
%     hold on;plot(c,d,'Color','k')
%     xlim([0 0.16])
%     ylim([0 0.16])
%     hold on;plot(mean(1./BestScoreOther),mean(1./BestScoreSelf),'+')
figure;hold on;
for jj=1:22
    plot(jj,BestSimX(:,jj),'*','Color','k')
    plot(jj,BestSimX(jj,jj),'*','Color','r')
end

% 5/22/2009
%%% YES!
% p=0.0118 2-sided - very significant
% Use ContingSim median method in order to diminish effects of outliers,
% which had a particularly strong effect on #3.
% BUT...median doesn't do well -- median(Score)
    for j=1:22
        NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
        for i=1:22
           NewTargs=Predict(i).onset+NormTargs;
           ind1=find(NewTargs>1);
           ind2=find(NewTargs(ind1)<1800);
           %CSs(i,j).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidAC,70)));
           Score2(i,j)=sum(abs(CSs(i,j).data(Predict(i).onset:Predict(i).onset+350)./max(CSs(i,j).data(Predict(i).onset:Predict(i).onset+350))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
            default=1./max(Prediction(i).data(Predict(i).onset:Predict(i).offset));
           for x=1:200
                fac(x)=default-0.8+0.005*x;
                matchscore(x)=sum(abs(fac(x)*CSs(i,j).data(Predict(i).onset:Predict(i).onset+350)./max(CSs(i,j).data(Predict(i).onset:Predict(i).onset+350))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
           end
           Score2(i,j)=min(matchscore);
            
        end
    end
    
    g=[1:28];
    for i=1:28
        indic1=find(Score(g~=i,i));
        indic2=find(Score(indic1,i)>0);
        BestScoreSelf2(i)=(Score(i,i));
        BestScoreOther2(i)=median(Score(indic1(indic2),i));
    end
    for i=1:22
        numX2(i)=length(find(BestSimX(1:22,i)<BestSimX(i,i)));
    end
    
     figure;hold on;
for i=1:25
    subplot(5,5,i);plot(ContingSim2(Predict(i).Targeting,Predict(i).ResidAC,80))
end



%%% 5/23/2009
%%%%%%%%%%%%%
%%%%%%%%%%%%%
% Parameter 1
%%%%%%%%%%%%%
%%%%%%%%%%%%%
for i=1:28
middle=Predict(i).onset+0.5*(Predict(i).offset-Predict(i).onset);
abc=(jccrossco2(Predict(i).ResidAC,middle));
hh(i,:)=abc(middle-150:middle+150)./max(abc(middle-150:middle+150));
end

for i=25
middle=Predict(i).onset+0.5*(Predict(i).offset-Predict(i).onset);
abc=(jccrossco2(Predict(i).ResidAC,middle));
hh(i,:)=abc(middle-150:middle+150)./max(abc(middle-150:middle+150));
end
% 118 is the center --- 150-32
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
% Parameter1 - precision of baseline variability
% Parameter2 - precision of targeting
% Parameter3 - precision of learning
% resids - precision of learning unaccounted for by precision of targeting

% NOTE THAT 70 and 170 are also significant AND avoid outliers.
parameter1=hh(:,70)+hh(:,170);
for i=1:28
parameter2(i)=std(Predict(i).Targeting);
end
for i=1:28
    [a,b]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    b=Predict(i).onset+b;
    parameter3(i)=Predict(i).LearnedNorm(b+80);
end
% for notes
for i=[4 8 9 12 16 17 21]
    [a,b]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
    b=Predict(i).onset+b;
    parameter3(i)=Predict(i).LearnedNorm(b-80);
end

p=polyfit(parameter2,parameter3,1);
f=polyval(p,parameter2);
resids=parameter3-f; % variation unexplained by targeting - negative values indicate 
                    % steeper than expected learning, positive values
                    % indicate broader than expected learning.
figure;plot(parameter1,resids,'*')
corrcoef(parameter1,resids)
%%%%%%%%%%%%
p2=polyfit(fwhm,parameter3,1);
f2=polyval(p2,fwhm);
resids2=parameter3-f2; % variation unexplained by targeting - negative values indicate 
                    % steeper than expected learning, positive values
                    % indicate broader than expected learning.
figure;plot(parameter2',resids2,'*')
corrcoef(parameter2',resids2)
% no relationship! - r=0.04

%%%%%%%%%%%% 
%%% More conventional way to get parameter1 -FWHM
inas=[1:28];
for ii=1:28
    k=inas(ii);
    clear ac
    middle=round(Predict(ii).onset+0.5*(Predict(ii).offset-Predict(ii).onset));
    for i=1:size(Predict(k).ResidAC,2)
        ac(i,:)=((xcorr(Predict(k).ResidAC(middle-80:middle+80,i))));
    end
    ccac(k).data=mean(ac)./max(mean(ac));
end
for ii=1:28
    indi=min(find(ccac(ii).data>0.5));
    fwhm(ii)=160-indi;
end
%%% This metric agrees - although FWat0.7 is much better than at 0.5
%%%%%%%%%%%%

%%%%% reinforced variability with this paradigm is generally consistent
%%%%% with structure of all variability
for i=1:28
    middle=Predict(i).onset+0.5*(Predict(i).offset-Predict(i).onset);
    abc=ContingSim2(middle,Predict(i).ResidAC,70);
    tt(i,:)=abc(middle-150:middle+150)./max(abc(middle-150:middle+150));
end
mmo=tt(:,70)+tt(:,230);
figure;plot(mmo(ind(1:18)),parameter3(ind(1:18)),'*')
%%% HIGHLY SIGNIFICANT!!!
% mmo vs. parameter3 resids - gives you best result and spreads out the
% data more - p<0.01!!!

% contingsim is more significant than crossco and thus generally is a better
% predictor - talk about how it includes the specific reinforced population
% - only contingsim does significantly better on the set w/o the outliers.

for i=1:17
selfers(i)=BestSims2(i,i);
medothers(i)=median(BestSims2(:,i));
mnothers(i)=mean(BestSims2(:,i));
end

for i=[1:18 20:22]
selfers(i)=Score(i,i);
medothers(i)=median(Score(1:22,i));
mnothers(i)=mean(Score(1:22,i));
end

%%% AC of INTACT vs LESION
inas=[3 7 10 13:22 24:28];
for ii=1:11
    k=inas(ii);
    clear ac
    clear ina
    for i=1:size(Predict(k).ResidAC,2)
        ac(i,:)=((xcorr(Predict(k).ResidAC(middle(k)-80:middle(k)+80,i))));
    end
    for i=1:size(Predict(k).ResidINA,2)
        ina(i,:)=((xcorr(Predict(k).ResidINA(middle(k)-80:middle(k)+80,i))));
    end
    ccac(k).data=mean(ac)./max(mean(ac));
    ccina(k).data=mean(ina)./max(mean(ina));
    dat(ii)=sum(ccac(k).data-ccina(k).data)/length(ccac(k).data);
end
for ii=indLes
    indiAC=min(find(ccac(ii).data>0.5));
    fwhmAC(ii)=160-indiAC;
    indiINA=min(find(ccina(ii).data>0.5));
    fwhmINA(ii)=160-indiINA;
end
for ii=[7 13 14 15 16 19 24:28]
    CVratio(ii)=1-mean(std(Predict(ii).ResidINA(Predict(ii).onset+50:Predict(ii).offset,:)'))./mean(std(Predict(ii).ResidAC(Predict(ii).onset+50:Predict(ii).offset,:)'));
end

for i=indLes
middle=Predict(i).onset+0.5*(Predict(i).offset-Predict(i).onset);
abc=(jccrossco2(Predict(i).ResidINA,middle));
hhLES(i,:)=abc(middle-150:middle+150)./max(abc(middle-150:middle+150));
end
parameter1LES=hhLES(:,70)+hhLES(:,170);

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
t=0:1:110;
indLes=[7 13 14 15 16 18 19 21 22 25 26 27 28];
%%% FIGURE 4 - LMAN influence %%%%%%%%%
figure;hold on;
% Fig.4A - example
subplot(3,3,1)
plot(Predict(25).ResidINA(240:430,[15 18 20 21 25]))
subplot(3,3,2)
plot(Predict(25).ResidAC(240:430,[15 18 20 21 25]))
subplot(3,3,3);hold on;
% LMAN INAs
plot(((fwhmINA(indLes([1:3 10:11]))*2)/8),((fwhmAC(indLes([1:3 10:11]))*2)/8),'*','Color','r')
% LMAN Lesions 
plot(((fwhmINA(indLes(4:9))*2)/8),((fwhmAC(indLes(4:9))*2)/8),'*','Color','b')
% RA AP-5
plot(((fwhmINA(indLes(12:13))*2)/8),((fwhmAC(indLes(12:13))*2)/8),'*','Color','g')
% How to make marker size give indication of SE???
plot(mean(fwhmINA(indLes)/4),mean(fwhmAC(indLes)/4),'+','MarkerSize',15,'Color','k')
plot(t,t,'-')
xlim([0 20]);ylim([0 20])

subplot(3,3,4)
Intact25=jccrossco(Predict(25).ResidAC,Predict(25).Targeting);
Les25=jccrossco(Predict(25).ResidINA,Predict(25).Targeting);
hold on;
plot(Les25(240:430)*1.7,'r');plot(Intact25(240:430)*1.15,'b');plot(Predict(25).LearnedNorm(240:430),'k')
subplot(3,3,5)
hold on;plot(IntactPred(indLes([1:3 10:11])),LesionPred(indLes([1:3 10:11])),'*','Color','r')
plot(IntactPred(indLes(4:9)),LesionPred(indLes(4:9)),'*','Color','b')
plot(IntactPred(indLes(12:13)),LesionPred(indLes(12:13)),'*','Color','g')
hold on;plot(t,t,'-','Color','k')
xlim([0 110]);ylim([0 110])
plot(mean(IntactPred(indLes)),mean(LesionPred(indLes)),'+','MarkerSize',15,'Color','k')



%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% This demonstrates how DMP has a fluctuating correlational
% structure, whereas DMP+AFP has a steadily decaying correlational
% structure.
    for i=indLes
    start=Predict(i).onset;
    abc=(jccrossco2(Predict(i).ResidINA,start+50));
    hh2LES(i,:)=abc(start+50:start+150)./max(abc(start+50:start+150));
    end
    for i=indLes
    start=Predict(i).onset;
    abc=(jccrossco2(Predict(i).ResidAC,start+50));
    hh2(i,:)=abc(start+50:start+150)./max(abc(start+50:start+150));
    end
    figure;plot(hh2LES','r')
    hold on;plot(hh2','k')

 %%%%
 %% Bimodal distribution
 % Width
 for i=1:28
     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
 end
 % Characterization of learning
         figure;hold on;
        abb=zeros(28,1400);
         for i=1:28
             [a,b]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
             dister1=b;
             dister2=notewidth(i)*8-b;

             abb(i,700-dister1:700)=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+b);
             abb(i,700:700+dister2)=Predict(i).LearnedNorm(Predict(i).onset+b:Predict(i).offset);
         end
        figure;plot(abb','*')
        %
        mnabb=zeros(1,1400);
        for i=1:1400
            ind=find(abb(:,i)>0);
            if ~isempty(ind)
            mnabb(i)=mean(abb(ind,i));
            end
        end
        hold on;plot(mnabb,'+')
  Figure1E.m
  
  %%%%%%%%%%%
  %%%%%%%%%%%
  % 6.04.09
 
  % FIGURE A1A1
  load FigA1.mat % located in /cardinal/FiguresA
        % createfigureA1A.m
        distance=112; % 30ms-16ms - distance between onset of spec and onset of contour
        initpt=250;
        finpt=900;
        % third harmonic
        meanBase=(3*mean(Alldatasecondseven1sig(7).baselineAC'));
        meanShifted=(3*0.5*(mean(Alldatasecondseven1sig(7).exp(17).selectedpitchcurves')+mean(Alldatasecondseven1sig(7).exp(18).selectedpitchcurves')));
        createfigureA1A(log(OR92heatmap.baseline))
        xlim([-0.01 0.083]);ylim([0 10000])
        % ADJUST limits on colormap to 4/15
        hold on;plot(t1(initpt-distance:770-distance),meanBase(initpt:770)/3,'LineWidth',2)
        %%
        createfigureA1A(log(OR92heatmap.shifted))
        xlim([-0.01 0.083]);ylim([0 10000])
        % ADJUST limits on colormap to 4/15
        hold on;plot(t1(initpt-distance:770-distance),meanShifted(initpt:770)/3,'LineWidth',2)
        %%
        createfigureA1A(log(OR92heatmap.baseline))
        xlim([-0.01 0.083]);ylim([6300 7800])
        hold on;plot(t1(initpt-distance:finpt-distance),meanBase(initpt:finpt),'LineWidth',2)
        %%
        createfigureA1A(log(OR92heatmap.shifted))
        xlim([-0.01 0.083]);ylim([6300 7800])
        % ADJUST limits on colormap to 4/15
        hold on;plot(t1(initpt-distance:finpt-distance),meanShifted(initpt:finpt),'LineWidth',2)
  % FIGURE A1A2
        prctshifted=(100*Predict(12).Learned)./mean(Alldatasecondseven1sig(7).baselineAC');
        figure;subplot(323);plot(prctshifted(initpt:finpt),'LineWidth',2,'Color','k')
        xlim([0 finpt-initpt])
        ylim([0 6])
                %   or92ht=[];
                %   for i=1:18 
                %      or92ht=[or92ht Alldatasecondseven1sig(7).exp(i).selectedpitchcurves];
                %   end
                %   for i=1:size(or92ht,2)
                %      or92htprct(:,i)=-100+(100*or92ht(:,i))./mean(Alldatasecondseven1sig(7).baselineAC')';
                %   end
                %   or92ht=or92ht(initpt:finpt,:);
          % 3 is a scaling factor that is arbitrary - get better
                %       or92htzscore=[];
                %       for i=1:size(or92ht,2)
                %           or92htzscore(:,i)=or92ht(:,i)'./std(Alldatasecondseven1sig(7).baselineAC(initpt:finpt,:)');
                %       end     
        subplot(3,2,1);image(or92htprct(initpt:finpt,:)'.*10);syn;
        % manually adjust y-axis - thousands of trials after 162(FB onset)
        [x1,y1]=hist(Predict(12).Targeting-initpt);
        subplot(325);stairs(y1,x1./length(Predict(12).Targeting),'LineWidth',2,'Color','k')
        xlim([0 finpt-initpt])
  % FIGURE A1B1
          figure;subplot(1,1,1,'XTickLabel',{'0','15','30','45'},'XTick',[0 120 240 360])
          hold on;
          plot(Predict(4).LearnedNorm(Predict(4).onset:Predict(4).offset),'LineWidth',2)
          plot(Predict(5).LearnedNorm(Predict(5).onset:Predict(5).offset),'LineWidth',2,'Color','k')
          plot(Predict(6).LearnedNorm(Predict(6).onset:Predict(6).offset),'LineWidth',2,'Color','r')
          [x1,y1]=hist(Predict(4).Targeting-Predict(4).onset,20);
          stairs(y1,x1./length(Predict(4).Targeting),'LineWidth',2)
          [x1,y1]=hist(Predict(5).Targeting-Predict(5).onset,20);
          stairs(y1,x1./length(Predict(5).Targeting),'LineWidth',2,'Color','k')
          [x1,y1]=hist(Predict(6).Targeting-Predict(6).onset,20);
          stairs(y1,x1./length(Predict(6).Targeting),'LineWidth',2,'Color','r')
          xlim([0 400])
          ylim([0 1.15])
          plot([median(Predict(4).Targeting)-Predict(4).onset median(Predict(4).Targeting)-Predict(4).onset],[0.1 0.3],'-')
          plot([median(Predict(5).Targeting)-Predict(5).onset median(Predict(5).Targeting)-Predict(5).onset],[0.1 0.3],'-','Color','k')
          plot([median(Predict(6).Targeting)-Predict(6).onset median(Predict(6).Targeting)-Predict(6).onset],[0.1 0.3],'-','Color','r')
          [a,b4]=max(Predict(4).LearnedNorm(Predict(4).onset:Predict(4).offset));
          [a,b5]=max(Predict(5).LearnedNorm(Predict(5).onset:Predict(5).offset));
          [a,b6]=max(Predict(6).LearnedNorm(Predict(6).onset:Predict(6).offset));
          plot([b4 b4],[0.9 1.1],'-')
          plot([b5 b5],[0.9 1.1],'-','Color','k')
          plot([b6 b6],[0.9 1.1],'-','Color','r')
  % FIGURE A1B2
        for i=1:28
           [a,tops(i)]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
           middle(i)=median(Predict(i).Targeting)-Predict(i).onset;
        end
        % get in milliseconds
        middle=middle/8;
        tops=tops/8;
        figure;
        plot(middle,tops,'.','MarkerSize',25,'Color','k')
        p=polyfit(middle,tops,1);
        t=0:0.1:100;
        hold on;plot(t,p(2)+p(1)*t,'b')
        xlim([0 70])
        ylim([0 70])
        
  %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%
  % TONS OF CODE BELOW -- this generates the global comparison of
  % variability compared with adaptation using variability similations of
  % 70th percentile and 50th percentile
  % Commented code centers on max adaptation
  % Uncommented (runnable) code centers on median targeting position
     % FIGURE A1C --- new A1E
                                                            %   
                                                            %                   % Mean variation across all experiments w/o targeting info
                                                            %                       figure;hold on;
                                                            %                       for i=1:28
                                                            %                           if isequal(Predict(i).direction,'up')
                                                            %                               prctile=70;
                                                            %                           else
                                                            %                               prctile=30;
                                                            %                           end
                                                            %                           CSraw(i).data=ContingSim2(median(Predict(i).Targeting),Predict(i).ResidAC,prctile);
                                                            %                           aax=CSraw(i).data;
                                                            %                           if prctile==70
                                                            %                               [a,btop]=max(aax(Predict(i).onset:Predict(i).offset));
                                                            %                           else
                                                            %                               [a,btop]=min(aax(Predict(i).onset:Predict(i).offset));
                                                            %                           end
                                                            %                           abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                            %                           left=length(abb);
                                                            %                           abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                            %                           right=length(abb)-left;
                                                            %                           t=-1*left:1:right-1;
                                                            %                           abb=abb/a;
                                                            %                           plot(t/8,abb,'Linewidth',2,'Color','k')
                                                            %                       end
                                                            %                          for i=1:28
                                                            %                              notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                            %                          end
                                                            %                          abb=zeros(28,1400);
                                                            %                          for i=1:28
                                                            % 
                                                            %                              [a,b]=max(abs(CSraw(i).data(Predict(i).onset:Predict(i).offset)));
                                                            %                              dister1=b;
                                                            %                              dister2=notewidth(i)*8-b;
                                                            % 
                                                            %                              abb(i,700-dister1:700)=abs(CSraw(i).data(Predict(i).onset:Predict(i).onset+b)/a);
                                                            %                              abb(i,700:700+dister2)=abs(CSraw(i).data(Predict(i).onset+b:Predict(i).offset)/a);
                                                            %                          end
                                                            % 
                                                            %                           mnabbNT=zeros(1,1400);
                                                            %                           seabbNT=zeros(1,1400);
                                                            %                           for i=1:1400
                                                            %                               ind=find(abb(:,i)>0);
                                                            %                               if ~isempty(ind)
                                                            %                                   mnabbNT(i)=mean(abb(ind,i));
                                                            %                                   seabbNT(i)=std(abb(ind,i))/sqrt(length(ind));
                                                            %                               end
                                                            %                           end
                                                            %                           t=-542:1:559;
                                                            %                       hold on;plot(t/8,mnabbNT(158:1259),'r','Linewidth',3)
                                                            % 
                                                            %                  % CS predictions - variability with targeting
                                                            %                       figure;hold on;
                                                            %                       for i=1:28
                                                            %                           aax=CSs2(i,i).data;
                        %                                                               if isequal(Predict(i).direction,'up')
                        %                                                                 [a,btop]=max(aax(Predict(i).onset:Predict(i).offset));
                        %                                                               else
                        %                                                                  [a,btop]=min(aax(Predict(i).onset:Predict(i).offset));
                        %                                                               end
                                                            %                           abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                            %                           left=length(abb);
                                                            %                           abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                            %                           right=length(abb)-left;
                                                            %                           t=-1*left:1:right-1;
                                                            %                           abb=abb/a;
                                                            %                           plot(t/8,abb,'Linewidth',2,'Color','k')
                                                            %                       end
                                                            %                          for i=1:28
                                                            %                              notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                            %                          end
                                                            %                          abb=zeros(28,1400);
                                                            %                          for i=1:28
                                                            %                              [a,b]=max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                                                            %                              dister1=b;
                                                            %                              dister2=notewidth(i)*8-b;
                                                            % 
                                                            %                              abb(i,700-dister1:700)=abs(CSs2(i,i).data(Predict(i).onset:Predict(i).onset+b)/a);
                                                            %                              abb(i,700:700+dister2)=abs(CSs2(i,i).data(Predict(i).onset+b:Predict(i).offset)/a);
                                                            %                          end
                                                            % 
                                                            %                           mnabbT=zeros(1,1400);
                                                            %                           seabbT=zeros(1,1400);
                                                            %                           for i=1:1400
                                                            %                               ind=find(abb(:,i)>0);
                                                            %                               if ~isempty(ind)
                                                            %                                   mnabbT(i)=mean(abb(ind,i));
                                                            %                                   seabbT(i)=std(abb(ind,i))/sqrt(length(ind));
                                                            %                               end
                                                            %                           end
                                                            %                           t=-542:1:559;
                                                            %                       hold on;plot(t/8,mnabbT(158:1259),'g','Linewidth',3)
                                                            % 
                                                            % 
                                                            %                   % Mean adaptation across all experiments
                                                            %                       figure;hold on
                                                            %                       for i=1:28
                                                            %                           [a,btop]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                                                            %                           abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                                                            %                           left=length(abb);
                                                            %                           abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                                                            %                           right=length(abb)-left;
                                                            %                           t=-1*left:1:right-1;
                                                            %                           plot(t/8,abb,'Linewidth',2,'Color','k')
                                                            %                       end
                                                            %                          for i=1:28
                                                            %                              notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                            %                          end
                                                            %                          abb=zeros(28,1400);
                                                            %                          for i=1:28
                                                            %                              [a,b]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                                                            %                              dister1=b;
                                                            %                              dister2=notewidth(i)*8-b;
                                                            % 
                                                            %                              abb(i,700-dister1:700)=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+b);
                                                            %                              abb(i,700:700+dister2)=Predict(i).LearnedNorm(Predict(i).onset+b:Predict(i).offset);
                                                            %                          end
                                                            % 
                                                            %                           mnabb=zeros(1,1400);
                                                            %                           seabb=zeros(1,1400);
                                                            %                           for i=1:1400
                                                            %                               ind=find(abb(:,i)>0);
                                                            %                               if ~isempty(ind)
                                                            %                                   mnabb(i)=mean(abb(ind,i));
                                                            %                                   seabb(i)=std(abb(ind,i))/sqrt(length(ind));
                                                            %                               end
                                                            %                           end
                                                            %                           t=-542:1:559;
                                                            %                       hold on;plot(t/8,mnabb(158:1259),'r','Linewidth',3)
                                                            %                       plot([-10 -10],[0 1.2],'-')
                                                            %                       plot([10 10],[0 1.2],'-')
                                                            %                       ylim([0 1.1])
                                                            %                       xlim([-70 70])


                                            % Mean adaptation across all experiments - centered at median targeting position
                                            % THIS GENERATES the New A1E plot
                                              figure;hold on
                                              for i=1:28
                                                  [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                                                  abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28
                                                     b=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                                                     b=round(b);
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+b);
                                                     abb(i,700:700+dister2)=Predict(i).LearnedNorm(Predict(i).onset+b:Predict(i).offset);
                                                 end
                                                  mnabb=zeros(1,1400);
                                                  seabb=zeros(1,1400);
                                                  for i=1:1400
                                                      ind=find(abb(:,i)>0);
                                                      if ~isempty(ind)
                                                          mnabb(i)=mean(abb(ind,i));
                                                          seabb(i)=std(abb(ind,i))/sqrt(length(ind));
                                                      end
                                                  end
                                                  t=-542:1:559;
                                              hold on;plot(t/8,mnabb(158:1259)/max(mnabb(158:1259)),'r','Linewidth',3)
                        %                       plot([-10 -10],[0 1.2],'-')
                        %                       plot([10 10],[0 1.2],'-')
                                              ylim([0 1.05])
                                              xlim([-40 40])
                                              
                               %%%%%%%%
                               %%%%%%%%%
                               %%% DISPROVE MODEL 1 - "Model1iswrong.eps"
                                            for i=1:28
                                                sd2targLL(i)=abb(i,700-round(std(Predict(i).Targeting))*2-80);
                                                sd2targL0(i)=abb(i,700-round(std(Predict(i).Targeting))*2);
                                                sd2targLR(i)=abb(i,700-round(std(Predict(i).Targeting))*2+80);
                                                sd2targRL(i)=abb(i,700+round(std(Predict(i).Targeting))*2-80);
                                                sd2targR0(i)=abb(i,700+round(std(Predict(i).Targeting))*2);
                                                sd2targRR(i)=abb(i,700+round(std(Predict(i).Targeting))*2+80);
                                            end
                                            figure;hold on
                                            count=0;
                                            for i=1:28
                                                if sd2targLL(i)>0
                                                    count=count+1;
                                                    plot(sd2targLR(i)-sd2targL0(i),sd2targL0(i)-sd2targLL(i),'.','MarkerSize',15,'Color','r')
                                                    xls(count)=sd2targLR(i)-sd2targL0(i);
                                                    yls(count)=sd2targL0(i)-sd2targLL(i);
                                                end
                                                
                                                if sd2targRR(i)>0
                                                    count=count+1;
                                                    plot(sd2targRL(i)-sd2targR0(i),sd2targR0(i)-sd2targRR(i),'.','MarkerSize',15)
                                                    xls(count)=sd2targRL(i)-sd2targR0(i);
                                                    yls(count)=sd2targR0(i)-sd2targRR(i);

                                                end
                                            end
                                            t=-1:0.1:1;
                                            hold on;plot(t,t,'Color','k')
                                            xlim([-0.15 0.35]);ylim([-0.15 0.35])

                                           % Mean variation across all experiments w/o targeting info - centered at targ position
                                              figure;hold on;
                                              for i=1:28
                                                  if isequal(Predict(i).direction,'up')
                                                      prctile=70;
                                                  else
                                                      prctile=30;
                                                  end
                                                  CSraw(i).data=ContingSim2(median(Predict(i).Targeting),Predict(i).ResidAC,prctile);
                                                  aax=CSraw(i).data;
                                                  if isequal(Predict(i).direction,'up')
                                                      a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                                                  else
                                                      a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                                                  end

                                                      btop=median(Predict(i).Targeting)-Predict(i).onset; %max(aax(Predict(i).onset:Predict(i).offset));

                                                  abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  abb=abb/a(i);
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28

                                                     b=median(Predict(i).Targeting)-Predict(i).onset; %max(abs(CSraw(i).data(Predict(i).onset:Predict(i).offset)));
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=abs(CSraw(i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                                                     abb(i,700:700+dister2)=abs(CSraw(i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                                                 end

                                                  mnabbNT=zeros(1,1400);
                                                  seabbNT=zeros(1,1400);
                                                  for i=1:1400
                                                      ind=find(abb(:,i)>0);
                                                      if ~isempty(ind)
                                                          mnabbNT(i)=mean(abb(ind,i));
                                                          seabbNT(i)=std(abb(ind,i))/sqrt(length(ind));
                                                      end
                                                  end
                                                  t=-542:1:559;
                                              hold on;plot(t/8,mnabbNT(158:1259),'r','Linewidth',3)

                                         % CS predictions - variability with targeting - centered at targ position 
                                              figure;hold on;
                                              for i=1:22
                                                  aax=CSs2(i,i).data;
                                                  if isequal(Predict(i).direction,'up')
                                                      a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                                                  else
                                                      a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                                                  end

                                                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                                                  abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  abb=abb/a(i);
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28
                                                     b=median(Predict(i).Targeting)-Predict(i).onset; %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=abs(CSs2(i,i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                                                     abb(i,700:700+dister2)=abs(CSs2(i,i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                                                 end

                                                  mnabbT=zeros(1,1400);
                                                  seabbT=zeros(1,1400);
                                                  for i=1:1400
                                                      ind=find(abb(:,i)>0);
                                                      if ~isempty(ind)
                                                          mnabbT(i)=mean(abb(ind,i));
                                                          seabbT(i)=std(abb(ind,i))/sqrt(length(ind));
                                                      end
                                                  end
                                                  t=-542:1:559;
                                              hold on;plot(t/8,mnabbT(158:1259),'g','Linewidth',3) 
                                                                  %%%%%
                                            %%% FINAL PLOT - generates "A1H"
                                            j1=max(mnabb(158:1259));
                                            j2=max(mnabbT(158:1000));
                                            j3=max(mnabbNT(158:1259));
                                             figure;hold on;plot(t/8,mnabb(158:1259)/j1,'b','LineWidth',2)
                                             plot([t/8;t/8],[mnabb(158:1259)/j1+seabb(158:1259)/j1;mnabb(158:1259)/j1-seabb(158:1259)/j1],'color','b')
                                            plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',3) % targ imprecision included
                                            plot(t/8,mnabbNT(158:1259)/j3,'k','LineWidth',3) % no targ imprecision
                                            xlim([-40 40]);ylim([0 1.05])
                                        %%%%%%%%
                                        %%%%%%%%%
                                        %%%%%%%%%%
                                        %%%%%%%
                        %%% Try 50/50 instead of 70/30
                                          % Mean variation across all experiments w/o targeting info - centered at targ position
                                              figure;hold on;
                                              for i=1:28
                                                  if isequal(Predict(i).direction,'up')
                                                      prctile=51;
                                                  else
                                                      prctile=49;
                                                  end
                                                  CSraw(i).data=ContingSim2(median(Predict(i).Targeting),Predict(i).ResidAC,prctile);
                                                  aax=CSraw(i).data;
                                                  if isequal(Predict(i).direction,'up')
                                                      a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                                                  else
                                                      a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                                                  end

                                                      btop=median(Predict(i).Targeting)-Predict(i).onset; %max(aax(Predict(i).onset:Predict(i).offset));

                                                  abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  abb=abb/a(i);
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28

                                                     b=median(Predict(i).Targeting)-Predict(i).onset; %max(abs(CSraw(i).data(Predict(i).onset:Predict(i).offset)));
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=abs(CSraw(i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                                                     abb(i,700:700+dister2)=abs(CSraw(i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                                                 end

                                                  mnabbNT=zeros(1,1400);
                                                  seabbNT=zeros(1,1400);
                                                  for i=1:1400
                                                      ind=find(abb(:,i)>0);
                                                      if ~isempty(ind)
                                                          mnabbNT(i)=mean(abb(ind,i));
                                                          seabbNT(i)=std(abb(ind,i))/sqrt(length(ind));
                                                      end
                                                  end
                                                  t=-542:1:559;
                                              hold on;plot(t/8,mnabbNT(158:1259),'r','Linewidth',3)

                                         % CS predictions - variability with targeting - centered at targ position 
                                              figure;hold on;
                                              for i=1:22
                                                  aax=CSs64(i,i).data;
                                                  if isequal(Predict(i).direction,'up')
                                                      a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                                                  else
                                                      a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                                                  end

                                                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                                                  abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  abb=abb/a(i);
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28
                                                     b=median(Predict(i).Targeting)-Predict(i).onset; %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=abs(CSs55(i,i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                                                     abb(i,700:700+dister2)=abs(CSs55(i,i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                                                 end

                                                  mnabbT=zeros(1,1400);
                                                  seabbT=zeros(1,1400);
                                                  for i=1:1400
                                                      ind=find(abb(:,i)>0);
                                                      if ~isempty(ind)
                                                          mnabbT(i)=mean(abb(ind,i));
                                                          seabbT(i)=std(abb(ind,i))/sqrt(length(ind));
                                                      end
                                                  end
                                                  t=-542:1:559;
                                              hold on;plot(t/8,mnabbT(158:1259),'g','Linewidth',3) 

                                            %%%%%
                                            %%% FINAL PLOT - generates "A1H"
                                            j1=max(mnabb(158:1259));
                                            j2=max(mnabbT(158:1000));
                                            j3=max(mnabbNT(158:1259));
                                             figure;hold on;plot(t/8,mnabb(158:1259)/j1,'b','LineWidth',2)
                                             plot([t/8;t/8],[mnabb(158:1259)/j1+seabb(158:1259)/j1;mnabb(158:1259)/j1-seabb(158:1259)/j1],'color','b')
                                            plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',3) % targ imprecision included
                                            plot(t/8,mnabbNT(158:1259)/j3,'k','LineWidth',3) % no targ imprecision
                                            xlim([-40 40]);ylim([0 1.05])

                                             %%% Plots for model description
                                            figure;plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',4);xlim([-40 40]);ylim([0 1.1])
                                             hold on;plot(t/8,[zeros(1,470) mnabbT(628:788)/j2 zeros(1,471)],'k','LineWidth',2)

                                             
                                             
                                       % A1F      
                                                       figure;hold on;
                                            subplot(241);plot(Predict(7).ResidAC(Predict(7).onset+30:Predict(7).onset+350,70:90))
                                                xlim([0 320]);ylim([-0.06 0.06]) 
                                            mtar=mean(Predict(7).Targeting)
                                            mprc=prctile(Predict(7).ResidAC(round(mtar),:),60)
                                            indH=find(Predict(7).ResidAC(round(mtar),:)>mprc)
                                            subplot(242);hold on;plot(Predict(7).ResidAC,'k')

                                            hold on;plot(Predict(7).ResidAC(:,indH),'r')
                                            hold on;plot(CSs64(7,7).data,'b','LineWidth',4)   
                                            xlim([Predict(7).onset+30 Predict(7).onset+350])
                                            ylim([-0.06 0.06])
                                            [xS,yS]=hist(Predict(7).Targeting-32,20);
                                            stairs(yS,-0.06+0.2*(xS/(2*length(Predict(7).Targeting))),'LineWidth',4,'Color','g')
                                           plot([0 1000],[0 0],'g','LineWidth',4)

                                            subplot(243);plot(Predict(7).LearnedNorm(Predict(7).onset+30:Predict(7).onset+350),'k','LineWidth',2)
                                                hold on;plot(abs(CSs64(7,7).data(Predict(7).onset+30:Predict(7).onset+350))/max(abs(CSs64(7,7).data(Predict(7).onset+30:Predict(7).onset+350)))*1.00,'b','LineWidth',2)

                                                t=0:1:350;
                                                xlim([0 320]);ylim([0 1.1])  

                                             
          %%%% SHOWS that w/in 40ms on either side, there is no bias toward
          %%%% broadening on the side with white noise
              seabb=zeros(1,1400);
              for i=1:1400
                  ind=find(abb(:,i)>0);
                  if ~isempty(ind)
                      seabb(i)=std(abb(ind,i))/sqrt(length(ind));
                  end
              end
              figure;hold on;plot(mnabb(700:1020)+seabb(700:1020),'LineWidth',2)
              plot(mnabb(700:1020)-seabb(700:1020),'LineWidth',2)
              t=320:-1:1;
              plot(t,mnabb(380:699)+seabb(380:699),'r','LineWidth',2)
              plot(t,mnabb(380:699)-seabb(380:699),'r','LineWidth',2)
              xlim([0 330])
              
         %%% Targ distn - Model 1
         
              
  % FIGURE A1D1
      % histogram of baseline song convolved with targeting distn
      for i=1:1000
          randpos1=round(rand*(length(Predict(12).Targeting)-1))+1;
          randpos2=round(rand*(size(Alldatasecondseven1sig(7).baselineAC,2)-1))+1;
          targdraw=Predict(12).Targeting(randpos1);
          pitchdraw(i)=mean(Alldatasecondseven1sig(7).baselineAC(targdraw-64:targdraw,randpos2));
      end
      [xS,yS]=hist(pitchdraw,20);
      figure;stairs(xS/1000,yS,'LineWidth',2)
      hold on;plot([0 0.3],[prctile(pitchdraw,70) prctile(pitchdraw,70)],'-','Color','k','LineWidth',2)
      xlim([0 0.2])
      figure;hist(pitchdraw,20) % change patch color to red
 
      % Get samples from folders
      %"avsamp" from 1119 - 1442.233
          figure;imagesc(OR92heatmap.tsamp,OR92heatmap.fsamp,log(OR92heatmap.avsamp));syn
            hold on;plot([1.857 1.865],[7700 7700],'-','Color','w','LineWidth',8)
            hold on;plot([0.061 0.069],[7700 7700],'-','Color','w','LineWidth',8)
            xlim([-0.4 2.4])
            ylim([0 10000])
            % ADJUST HEAT MAP to 'hot' and low/high to 4/15 - this yields
            % perfection!!!!
    %%% FIGURE A2
        % Long example - #15
        % Short example - #28
        % See lines ~1200-1300
        figure;hold on;
        % spectrograms
        subplot(3,5,1);imagesc(Predict(15).t,Predict(15).f,log(Predict(15).avnote));syn
        ylim([0 10000]); xlim([0.04 0.17])
        subplot(3,5,6);imagesc(Predict(28).t,Predict(28).f,log(Predict(28).avnote));syn
        ylim([0 10000]); xlim([-0.03 0.1])
        % Bimodal distribution
                for i=1:28
                    mlength(i)=Predict(i).offset-Predict(i).onset;
                end
        
                    %                 Alldata1sig=[Alldatafirstten1sig Alldatasecondseven1sig];
                    %                 for i=1:17
                    %                 mmean(i)=mean(mean(Alldata1sig(i).baselineAC(Alldata1sig(i).startnote:Alldata1sig(i).endnote,:)));
                    %                 end
                    %                 mmean1(1)=mmean(1);
                    %                 mmean1(2)=mmean(3);
                    %                 mmean1(23)=mmean(4);
                    %                 mmean1(24)=mmean(5);
                    %                 mmean1(3)=mmean(6);
                    %                 mmean1(4:6)=mmean(7:9);
                    %                 mmean1(7:12)=mmean(12:17);
                    %                 mmean1(25:26)=mmean(10:11);
                    %                 mmean=mmean1;
                    %                 mmean(13)=mean(mean(meancurves.bk61w42(Predict(13).onset:Predict(13).offset,:)));
                    %                 mmean(14)=mean(mean(meancurves.bk63w43(Predict(14).onset:Predict(14).offset,:)));
                    %                 mmean(16)=mean(mean(meancurves.pk42(Predict(16).onset:Predict(16).offset,:)));
                    %                 mmean(18)=mean(mean(meancurves.pk28(Predict(18).onset:Predict(18).offset,:)));
                    %                 mmean(19)=mean(mean(meancurves.bk13(Predict(19).onset:Predict(19).offset,:)));
                    %                 mmean(21)=mean(mean(meancurves.o50(Predict(21).onset:Predict(21).offset,:)));
                    %                 mmean(15)=mean(mean(meancurves.pk30(Predict(15).onset:Predict(15).offset,:)));
                    %                 mmean(17)=mean(mean(meancurves.pk29(Predict(17).onset:Predict(17).offset,:)));
                    %                 mmean(20)=mean(mean(meancurves.bk72(Predict(20).onset:Predict(20).offset,:)));
                    %                 mmean(22)=mean(mean(meancurves.o11(Predict(22).onset:Predict(22).offset,:)));
                    %                 mmean(27)=mean(mean(meancurves.pu57(Predict(27).onset:Predict(27).offset,:)));
                    %                 mmean(28)=mean(mean(meancurves.pu56(Predict(28).onset:Predict(28).offset,:)));
        subplot(3,5,11);plot(mlength,mmean,'.','MarkerSize',15)
        hold on;plot(mlength(23:28),mmean(23:28),'.','MarkerSize',15,'Color','r')
        xlim([80 680]);ylim([2000 4000])

        % A2C - learning and targeting - centered around peak
        tx=74:1:294;
        subplot(3,5,3);plot(tx,Predict(28).LearnedNorm(Predict(28).onset:Predict(28).offset),'r','LineWidth',2)
        hold on;plot(Predict(15).LearnedNorm(Predict(15).onset:Predict(15).offset),'LineWidth',2)
        [xS,yS]=hist(Predict(28).Targeting-Predict(28).onset+74,20);
        hold on;stairs(yS,xS/length(Predict(28).Targeting),'LineWidth',2,'Color','r')
        [xS,yS]=hist(Predict(15).Targeting-Predict(15).onset,20);
        hold on;stairs(yS,xS/length(Predict(15).Targeting),'LineWidth',2,'Color','b')
        plot([90 90],[0.5 1],'-','Color','k')
        plot([250 250],[0.5 1],'-','Color','k')
        ylim([0 1.1])
        xlim([0 400])

        % NOTE THAT I'm only removing the residuals with ugly transients
            % A2D - 20 residuals for long note - middle 25ms
                subplot(3,5,2);plot(Predict(15).ResidAC(800:1000,[5:15 17:25]))
                xlim([0 200]);ylim([-0.06 0.06])
            % A2E - 20 residuals for short note - middle 25ms
                subplot(3,5,7);plot(Predict(28).ResidAC(200:400,[36:37 39:41 45:59]))
                xlim([0 200]);ylim([-0.06 0.06])
        
        % A2E - autocorrelation metric
            subplot(3,5,8);plot(ccac(15).data,'LineWidth',2)
            hold on;plot(ccac(28).data,'r','LineWidth',2)
            plot([0 320],[0.5 0.5],'-','Color','k')
            xlim([0 320])
            ylim([0 1.1])
        % A2F - imprecision of targeting std (ms) vs. imprecision/generalization of adaptation (remaining after 10ms)
            subplot(3,5,5);plot(parameter2(1:22)/8,parameter3(1:22),'.','MarkerSize',15,'Color','b')
            hold on;plot(parameter2(23:28)/8,parameter3(23:28),'.','MarkerSize',15,'Color','r')
            p=polyfit(parameter2/8,parameter3,1);
            f=polyval(p,parameter2/8);
            resids1=parameter3-f; % variation unexplained by targeting - negative values indicate 
            t=0:1:50
            hold on;plot(t,p(2)+t*p(1),'Color','k')
            xlim([0 30])
            % r=0.25
            % m=0.0035
            
        % A2G - imprecision of variability (fwhm of autocorr in ms) vs. imprecision/generalization of adaptation (remaining after 10ms)
            subplot(3,5,4);plot(fwhm(1:22)/4,parameter3(1:22),'.','MarkerSize',15,'Color','b')
            hold on;plot(fwhm(23:28)/4,parameter3(23:28),'.','MarkerSize',15,'Color','r')
            p=polyfit(fwhm/4,parameter3,1);
            f=polyval(p,fwhm/4);
            resids2=parameter3-f; % variation unexplained by targeting - negative values indicate 
            t=0:1:50
            hold on;plot(t,p(2)+t*p(1),'Color','k')
            xlim([9 18])
            % r=0.42
            % m=0.017
           
        % A2H - targeting vs. resids from autocorr
            subplot(3,5,10);plot(parameter2(1:22)/8,resids2(1:22),'.','MarkerSize',15,'Color','b')
            hold on;plot(parameter2(23:28)/8,resids2(23:28),'.','MarkerSize',15,'Color','r')
            p=polyfit(parameter2/8,resids2,1);
            f=polyval(p,parameter2/8);
            t=0:1:50
            hold on;plot(t,p(2)+t*p(1),'Color','k')
            xlim([0 30])
            % r=0.10
            % m=-0.0123
            
        % A2I - autocorr vs resids from targeting
            subplot(3,5,9);plot(fwhm(1:22)/4,resids1(1:22),'.','MarkerSize',15,'Color','b')
            hold on;plot(fwhm(23:28)/4,resids1(23:28),'.','MarkerSize',15,'Color','r')
            p=polyfit(fwhm/4,resids1,1);
            f=polyval(p,fwhm/4);
            t=0:1:50
            hold on;plot(t,p(2)+t*p(1),'Color','k')
            xlim([9 18])
            % r=0.33
            % m=0.0133
            
            
    %%% FIGURE A3
            % Make the CS predictions for long notes 
                for j=1:22
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    if isequal('up',Predict(j).direction);value=70;else value=30;end
                    for i=1:22
                       NewTargs=Predict(i).onset+NormTargs;
                       ind1=find(NewTargs>1);
                       ind2=find(NewTargs(ind1)<1800);
                       CSs2(i,j).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidAC,value)));
                    end
                end
            % Make the CS predictions for short notes 
                for j=23:28
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    if isequal('up',Predict(j).direction);value=70;else value=30;end
                    for i=23:28
                        NewTargs=Predict(i).onset+NormTargs;
                        ind1=find(NewTargs>1);
                        ind2=find(NewTargs(ind1)<1800);
                        CSs2(i,j).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidAC,value)));
                    end
                end
            % Fit the CS predictions for long notes - optimal scaling factor
                for i=1:22
                    for j=1:22
                        count=0;
                        for k=0.5:0.01:2
                            count=count+1;
                            choices(count)=sum(abs((k*abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                        end
                        disters2(i,j)=min(choices);
                    end
                end
            % Fit the CS predictions for short notes - optimal scaling factor
                for i=23:28
                    for j=23:28
                        count=0;
                        for k=0.5:0.01:2
                            count=count+1;
                            choices(count)=sum(abs((k*abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+160))./max(abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                        end
                        disters2(i,j)=min(choices);
                    end
                end           
            % Plot the *unscaled* CS predictions against the actual learning
                for i=1:28
                    subplot(5,6,i);hold on;plot(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).onset+350))))
                    plot(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350),'r')
                    xlim([0 350]);ylim([0 1.1])
                end
            % DC - how good is this prediction
                    % DC - goodness of fit 
                    for i=1:22
                        BestScoreDC(i)=sum(abs(mean(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350))-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350)));
                    end
                    for i=23:28
                        BestScoreDC(i)=sum(abs(mean(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+160))-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+160)));
                    end                 
            % Self/Other - how good are these predictions
                for i=1:28
                    BestSelfCS(i)=disters2(i,i);
                end
                % long - other
                g=[1:22];
                for i=1:22
                    indn=find(g~=i);
                    BestOtherCS(i)=mean(disters2(indn,i));
                end
                % short - other
                g=[23:28];
                for i=23:28
                    indn=find(g~=i);
                    BestOtherCS(i)=mean(disters2(g(indn),i));
                end
                
                % normalize scores
                    BestSelfCSnorm(1:22)=(BestSelfCS(1:22)/350);
                    BestSelfCSnorm(23:28)=(BestSelfCS(23:28)/160);
                    BestOtherCSnorm(1:22)=(BestOtherCS(1:22)/350);
                    BestOtherCSnorm(23:28)=(BestOtherCS(23:28)/160);
                    BestScoreDCnorm(1:22)=(BestScoreDC(1:22)/350);
                    BestScoreDCnorm(23:28)=(BestScoreDC(23:28)/160);
%%% 7.19.09 %%% SAME SYLLABLE VS. OTHER SYLLABLES
    % exps 1and2, exps 4,5,6, exps 8and9
    [h,p]=ttest([mean(BestSelfCSnorm(1:2)) BestSelfCSnorm(3) mean(BestSelfCSnorm(4:6)) BestSelfCSnorm(7) mean(BestSelfCSnorm(8:9)) BestSelfCSnorm(10:28)],[mean(BestOtherCSnorm(1:2)) BestOtherCSnorm(3) mean(BestOtherCSnorm(4:6)) BestOtherCSnorm(7) mean(BestSelfCSnorm(8:9)) BestOtherCSnorm(10:28)])
    p=0.0227; %thus significant even for 70/30
%%%% Try 40/60
                    % Fit the CS predictions for long notes - optimal scaling factor
                        for i=1:22
                            for j=1:22
                                count=0;
                                for k=0.5:0.01:2
                                    count=count+1;
                                    choices(count)=sum(abs((k*abs(CSs64(i,j).data(Predict(i).onset+30:Predict(i).onset+350))./max(abs(CSs64(i,j).data(Predict(i).onset+30:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset+30:Predict(j).onset+350)));
                                end
                                disters2(i,j)=min(choices);
                            end
                        end
                    % Fit the CS predictions for short notes - optimal scaling factor
                        for i=23:28
                            for j=23:28
                                count=0;
                                for k=0.5:0.01:2
                                    count=count+1;
                                    choices(count)=sum(abs((k*abs(CSs64(i,j).data(Predict(i).onset:Predict(i).onset+160))./max(abs(CSs64(i,j).data(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                                end
                                disters2(i,j)=min(choices);
                            end
                        end           
                    % Plot the *unscaled* CS predictions against the actual learning
                        for i=1:28
                            subplot(5,6,i);hold on;plot(abs(CSs64(i,i).data(Predict(i).onset+30:Predict(i).onset+350))./max(abs(CSs64(i,i).data(Predict(i).onset:Predict(i).onset+350))))
                            plot(Predict(i).LearnedNorm(Predict(i).onset+30:Predict(i).onset+350),'r')
                            xlim([0 350]);ylim([0 1.1])
                        end
                    % DC - how good is this prediction
                            % DC - goodness of fit 
                            for i=1:22
                                BestScoreDC(i)=sum(abs(mean(Predict(i).LearnedNorm(Predict(i).onset+30:Predict(i).onset+350))-Predict(i).LearnedNorm(Predict(i).onset+30:Predict(i).onset+350)));
                            end
                            for i=23:28
                                BestScoreDC(i)=sum(abs(mean(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+160))-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+160)));
                            end                 
                    % Self/Other - how good are these predictions
                        for i=1:28
                            BestSelfCS(i)=disters2(i,i);
                        end
                        % long - other
                        g=[1:22];
                        for i=1:22
                            indn=find(g~=i);
                            BestOtherCS(i)=mean(disters2(indn,i));
                        end
                        % short - other
                        g=[23:28];
                        for i=23:28
                            indn=find(g~=i);
                            BestOtherCS(i)=mean(disters2(g(indn),i));
                        end

                        % normalize scores
                            BestSelfCSnorm64(1:22)=(BestSelfCS(1:22)/320);
                            BestSelfCSnorm64(23:28)=(BestSelfCS(23:28)/160);
                            BestOtherCSnorm64(1:22)=(BestOtherCS(1:22)/320);
                            BestOtherCSnorm64(23:28)=(BestOtherCS(23:28)/160);
                            BestScoreDCnorm64(1:22)=(BestScoreDC(1:22)/320);
                            BestScoreDCnorm64(23:28)=(BestScoreDC(23:28)/160);
    
                           [h,p]=ttest(BestSelfCSnorm64,BestOtherCSnorm64) 
                            [h,p]=ttest([mean(BestSelfCSnorm64(1:2)) BestSelfCSnorm64(3) mean(BestSelfCSnorm64(4:6)) BestSelfCSnorm64(7) mean(BestSelfCSnorm64(8:9)) BestSelfCSnorm64(10:28)],[mean(BestOtherCSnorm64(1:2)) BestOtherCSnorm64(3) mean(BestOtherCSnorm64(4:6)) BestOtherCSnorm64(7) mean(BestSelfCSnorm64(8:9)) BestOtherCSnorm64(10:28)])
                           [h,p]=ttest(BestSelfCSnorm64,BestScoreDCnorm64)
          % FIGURE - plot examples and statistics - FINAL!!!
            figure;hold on;
            subplot(241);plot(Predict(1).ResidAC(Predict(1).onset+30:Predict(1).onset+350,40:60))
                xlim([0 320]);ylim([-0.06 0.06]) 
            subplot(245);plot(Predict(4).ResidAC(Predict(4).onset+30:Predict(4).onset+350,70:90))
                xlim([0 320]);ylim([-0.06 0.06]) 
            mtar=mean(Predict(4).Targeting)
            mprc=prctile(Predict(4).ResidAC(round(mtar),:),60)
            indH=find(Predict(4).ResidAC(round(mtar),:)>mprc)
            subplot(246);plot(Predict(4).ResidAC,'k')
            hold on;plot(Predict(4).ResidAC(:,indH),'r')
            hold on;plot(CSs64(4,4).data,'b','LineWidth',4)   
            xlim([Predict(4).onset+30 Predict(4).onset+350])
            ylim([-0.06 0.06])
            [xS,yS]=hist(Predict(4).Targeting-32,20);
            stairs(yS,-0.06+0.2*(xS/length(Predict(4).Targeting)),'LineWidth',4,'Color','g')

            mtar=mean(Predict(1).Targeting)
            mprc=prctile(Predict(1).ResidAC(round(mtar),:),60)
            indH=find(Predict(1).ResidAC(round(mtar),:)>mprc)
            subplot(242);plot(Predict(1).ResidAC,'k')
            hold on;plot(Predict(1).ResidAC(:,indH),'r')
            hold on;plot(CSs64(1,1).data,'b','LineWidth',4)   
            xlim([Predict(1).onset+30 Predict(1).onset+350])
            ylim([-0.06 0.06])
            [xS,yS]=hist(Predict(1).Targeting-32,20);
            stairs(yS,-0.06+0.1*(xS/length(Predict(1).Targeting)),'LineWidth',4,'Color','g')
%%%% Note - if i do 70/30 instead - substitute exp #7 for exp#1
            subplot(243);plot(Predict(1).LearnedNorm(Predict(1).onset+30:Predict(1).onset+350),'k','LineWidth',2)
                hold on;plot(abs(CSs64(1,1).data(Predict(1).onset+30:Predict(1).onset+350))/max(abs(CSs64(1,1).data(Predict(1).onset+30:Predict(1).onset+350)))*1.00,'b','LineWidth',2)
                plot(abs(CSs64(4,1).data(Predict(4).onset+30:Predict(4).onset+350))/max(abs(CSs64(4,1).data(Predict(4).onset+30:Predict(4).onset+350)))*1.01,'r','LineWidth',2)   
                t=0:1:350;
                plot(t,mean(Predict(1).LearnedNorm(Predict(1).onset+30:Predict(1).onset+350)),'g','LineWidth',2)
                xlim([0 320]);ylim([0 1.1])
            subplot(247);plot(Predict(4).LearnedNorm(Predict(4).onset+30:Predict(4).onset+350),'k','LineWidth',2)
                hold on;plot(abs(CSs64(4,4).data(Predict(4).onset+30:Predict(4).onset+350))/max(abs(CSs64(4,4).data(Predict(4).onset+30:Predict(4).onset+350)))*1.05,'b','LineWidth',2)
                plot(abs(CSs64(1,4).data(Predict(1).onset+30:Predict(1).onset+350))/max(abs(CSs64(1,4).data(Predict(1).onset+30:Predict(1).onset+350)))*0.99,'r','LineWidth',2)   
                t=0:1:350;
                plot(t,mean(Predict(4).LearnedNorm(Predict(4).onset:Predict(4).onset+350)),'g','LineWidth',2)
                xlim([0 320]);ylim([0 1.1]) 
            subplot(248);hold on;
            plot(2,BestSelfCSnorm64,'.','MarkerSize',20,'Color','b')
                plot(2,mean(BestSelfCSnorm64),'.','MarkerSize',50,'Color','k')
                hold on;plot([1 2],[BestOtherCSnorm64;BestSelfCSnorm64],'-','Color','k')
                hold on;plot(1,[BestOtherCSnorm64],'.','MarkerSize',20,'Color','r')
                plot(1,mean(BestOtherCSnorm64),'.','MarkerSize',50,'Color','k')
                hold on;plot(3,[BestScoreDCnorm64],'.','MarkerSize',20,'Color','g')
                plot(3,mean(BestScoreDCnorm64),'.','MarkerSize',50,'Color','k')
                hold on;plot([2 3],[BestSelfCSnorm64;BestScoreDCnorm64],'-','Color','k')
                xlim([0.7 3.3]);ylim([0 0.25])
                    % coefficients of best fit for the plots above determined by the for loop below
                    % 7,7=1.03 % 4,4=1.05 % 7,4=0.94 % 4,7=1.01
                            for i=4
                                for j=7
                                    count=0;
                                    a=0.5:0.01:2;
                                    for k=a
                                        count=count+1;
                                        choices(count)=sum(abs((k*abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                                    end
                                    [b,c]=min(choices);
                                    a(c)
                                end
                            end
          % Control that's good to know about but not all that helpful:  
             % Wrong Direction generally doesn't matter except in 3 instances (#3,7,11)
                for j=1:22
                    if isequal('up',Predict(j).direction)
                        value=30;
                    else
                        value=70;
                    end
                    CSwrongdir(j).data=(mean(ContingSim2(Predict(j).Targeting,Predict(j).ResidAC,value)));
                end
                for i=1:22
                    count=0;
                    for k=0.5:0.01:2
                        count=count+1;
                        choices(count)=sum(abs((k*abs(CSwrongdir(i).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSwrongdir(i).data(Predict(i).onset:Predict(i).offset))))-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350)));
                    end
                    distersWrongDir(i)=min(choices);
                end
            %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    % CROSSCO - Here is a good figure for crossco method
                    %     prediction11_7=jccrossco(Predict(11).ResidAC,700+Predict(7).Targeting-Predict(7).onset);
                    %     prediction7_11=jccrossco(Predict(7).ResidAC,340+Predict(11).Targeting-Predict(11).onset);

                    % These two were chosen because both are targeted in the same region with
                    % precise targeting distributions, because differences in the structure of
                    % variability can be followed by eye to differences in prediction and
                    % predict subtle but apparent differences in learning.
                    % prediction1_15=jccrossco(Predict(1).ResidAC,Predict(1).onset+Predict(15).Targeting-Predict(15).onset);
                    % prediction15_1=jccrossco(Predict(15).ResidAC,Predict(15).onset+Predict(1).Targeting-Predict(1).onset);

                    figure;hold on;subplot(232);hold on
                    plot(Predict(15).LearnedNorm(Predict(15).onset:Predict(15).onset+350),'k','LineWidth',2)
                    plot(Predict(15).crosscoAC(Predict(15).onset:Predict(15).onset+350)/max(Predict(15).crosscoAC(Predict(15).onset:Predict(15).onset+350)),'b','LineWidth',2)
                    plot(prediction1_15(Predict(1).onset:Predict(1).onset+350)/max(prediction1_15(Predict(1).onset:Predict(1).onset+350)),'r','LineWidth',2)
                    t=0:1:350;
                    plot(t,mean(Predict(15).LearnedNorm(Predict(15).onset:Predict(15).onset+350)),'g','LineWidth',2)
                    xlim([0 350]); ylim([0 1.1])
                    subplot(235);hold on
                    plot(Predict(1).LearnedNorm(Predict(1).onset:Predict(1).onset+350),'k','LineWidth',2)
                    plot(Predict(1).crosscoAC(Predict(1).onset:Predict(1).onset+350)/max(Predict(1).crosscoAC(Predict(1).onset:Predict(1).onset+350)),'b','LineWidth',2)
                    plot(prediction15_1(Predict(15).onset:Predict(15).onset+350)/max(prediction15_1(Predict(15).onset:Predict(15).onset+350)),'r','LineWidth',2)
                    plot(t,mean(Predict(1).LearnedNorm(Predict(1).onset:Predict(1).onset+350)),'g','LineWidth',2)
                    xlim([0 350]); ylim([0 1.1])

                    subplot(231);plot(Predict(15).ResidAC(Predict(15).onset:Predict(15).onset+350,41:60))
                    xlim([0 350]);ylim([-0.06 0.06])
                    subplot(234);plot(Predict(1).ResidAC(Predict(1).onset:Predict(1).onset+350,[40:44 50:53 60:65 90:94]))
                    xlim([0 350]);ylim([-0.06 0.06])

                    %             figure;hold on;plot(1,1./(BestScoreNOTARGSelf/350),'*','Color','k')
                    %             plot(1,median(1./(BestScoreNOTARGSelf/350)),'.','MarkerSize',50,'Color','r')
                                subplot(236);hold on;plot(2,1./(BestScoreSelf/350),'.','MarkerSize',20,'Color','b')
                                plot(2,median(1./(BestScoreSelf/350)),'.','MarkerSize',50,'Color','k')
                    %             hold on;plot([1 2],[1./(BestScoreNOTARGSelf/350);1./(BestScoreSelf/350)],'-','Color','k')
                                hold on;plot([2 3],[1./(BestScoreSelf/350);1./(BestScoreOther/350)],'-','Color','k')
                                hold on;plot(3,[1./(BestScoreOther/350)],'.','MarkerSize',20,'Color','r')
                                plot(3,median(1./(BestScoreOther/350)),'.','MarkerSize',50,'Color','k')
                                hold on;plot(1,[1./(BestScoreDC/350)],'.','MarkerSize',20,'Color','g')
                                plot(1,median(1./(BestScoreDC/350)),'.','MarkerSize',50,'Color','k')
                                hold on;plot([1 2],[1./(BestScoreDC/350);1./(BestScoreSelf/350)],'-','Color','k')
                                xlim([0.7 3.3])

    % CROSSCO - Determine goodness of prediction  --- crossco method              
                    for j=1:22
                    %NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                        for i=1:22
                            NewTargs=Predict(i).onset+NormTargs;
                            %clear crossco;
                            %Prediction(i,j).data=jccrossco(Predict(i).ResidAC,NewTargs(find(NewTargs>64))); % round?
                            % best prediction
                            default=1./max(Prediction(i,j).data(Predict(i).onset:Predict(i).offset));
                            for x=1:200
                                fac(x)=default-0.8+0.005*x;
                                matchscore(x,:)=sum(abs(Prediction(i,j).data(Predict(i).onset:Predict(i).onset+350)*fac(x)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                            end
                            [BestSimX(i,j),ind]=min(matchscore);
                        end
                    end
            
                    BestSimX(3,20)=mean(BestSimX([1:2 4:22],20)); % NaN replaced
                    figure;plot(BestScoreSelf/350)
                    hold on;plot(BestScoreOther/350,'r')
                    g=[1:22];
                    for i=1:22
                        BestScoreSelf(i)=mean(BestSimX(i,i));
                        BestScoreOther(i)=mean(BestSimX(find(g~=i),i));
                    end
     % CROSSCO - Goodness of prediction w/o targeting
                    for j=1:22
                        NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                        for i=1:22
                            NewTargs=Predict(i).onset+NormTargs;
                            %clear crossco;
                            PD=jccrossco2(Predict(i).ResidAC,median(NewTargs)); % round?
                            % best prediction
                            default=1./max(Prediction(i,j).data(Predict(i).onset:Predict(i).offset));
                            for x=1:200
                                fac(x)=default-0.8+0.005*x;
                                matchscore(x,:)=sum(abs(PD(Predict(i).onset:Predict(i).onset+350)*fac(x)-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                            end
                            [BestSimNOTARG(i,j),ind]=min(matchscore);
                        end
                    end
                    BestSimNOTARG(3,20)=mean(BestSimNOTARG([1:2 4:22],20)); % NaN replaced
                    figure;plot(BestScoreSelf/350)
                    hold on;plot(BestScoreOther/350,'r')
                    g=[1:22];
                    for i=1:22
                        BestScoreNOTARGSelf(i)=mean(BestSimNOTARG(i,i));
                        BestScoreNOTARGOther(i)=mean(BestSimNOTARG(find(g~=i),i));
                    end
            

    % CROSSCO and CS methods - old examples
            
                                    % A3Ai - experiment 1, 20 residuals
                %             figure;plot(Predict(1).ResidAC(218:568,[20:31 161:170]))
                        % A3Ai - experiment 15, 20 residuals 
                        % Argument is that there is an obvious difference between
                        % these 3 patterns of variability - should be diff in predictions?
                        % We predict the prediction from 15 will be better for 15.
                            figure;subplot(322);plot(Predict(15).ResidAC(750:1100,1:20),'b')
                            subplot(321);plot(Predict(19).ResidAC(300:650,21:40),'r')
                            %figure;plot(Predict(7).ResidAC(340:690,61:80))

                            medtar=median(Predict(15).Targeting)-Predict(15).onset; % 145
                % %           Illustrative predictions based on just one point
                % %             prediction15_7a=jccrossco2(Predict(7).ResidAC,340+round(medtar));
                % %             prediction15_7b=jccrossco(Predict(7).ResidAC,340+Predict(15).Targeting-Predict(15).onset);
                %                 prediction15_15a=jccrossco2(Predict(15).ResidAC,750+round(medtar));             
                %                 prediction15_15b=Predict(15).crosscoAC;
                %                 prediction15_19a=jccrossco2(Predict(19).ResidAC,300+round(medtar));
                %                 prediction15_19b=jccrossco(Predict(19).ResidAC,300+Predict(15).Targeting-Predict(15).onset);
                %                 prediction19_19b=Predict(19).crosscoAC;
                %                 prediction19_15b=jccrossco(Predict(15).ResidAC,750+Predict(19).Targeting-Predict(19).onset);
                            subplot(323);hold on
                                plot(Predict(19).LearnedNorm(300:650),'k','LineWidth',2)
                                plot(prediction19_19b(300:650)./max(prediction19_19b(300:650)),'r','LineWidth',2)
                                plot(prediction19_15b(750:1100)./max(prediction19_15b(750:1100)),'b','LineWidth',2)
                                xlim([0 350]);ylim([0 1.1])
                            subplot(324);hold on
                                plot(Predict(15).LearnedNorm(750:1100),'k','LineWidth',2)
                                plot(prediction15_15b(750:1100)./max(prediction15_15b(750:1100)),'LineWidth',2)
                                plot(prediction15_19b(300:650)./max(prediction15_19b(300:650)),'r','LineWidth',2)
                                xlim([0 350]);ylim([0 1.1])
                %             %abc7_15=ContingSim2(Predict(15).Targeting-410-32,Predict(7).ResidAC,70);
                %                 abc19_15=ContingSim2(Predict(15).Targeting-450-32,Predict(19).ResidAC,70);
                %                 abc19_15=mean(abc19_15);
                %                 abc15=ContingSim2(Predict(15).Targeting-32,Predict(15).ResidAC,70);
                %                 abc15=mean(abc15);
                %                 abc19=ContingSim2(Predict(19).Targeting-32,Predict(19).ResidAC,70);
                %                 abc19=mean(abc19);
                %                 abc15_19=ContingSim2(Predict(19).Targeting+450-32,Predict(15).ResidAC,70);
                %                 abc15_19=mean(abc15_19);
                            subplot(325);hold on
                                plot(Predict(19).LearnedNorm(300:650),'k','LineWidth',2)
                                plot(abc15_19(750:1100)./max(abc15_19(750:1100)),'b','LineWidth',2)
                                plot(abc19(300:650)./max(abc19(300:650)),'r','LineWidth',2)
                                xlim([0 350]);ylim([0 1.1])
                            subplot(326);hold on;
                                plot(Predict(15).LearnedNorm(750:1100),'k','LineWidth',2)
                                plot(abc19_15(300:650)./max(abc19_15(300:650)),'r','LineWidth',2)
                                plot(abc15(750:1100)./max(abc15(750:1100)),'LineWidth',2)
                                xlim([0 350]);ylim([0 1.1])

% Here is an old example from crossco method.
                % %     figure;hold on;subplot(222)
                % %         plot(prediction7_11(340:690)/max(prediction7_11(340:690)),'r','LineWidth',2)
                % %         hold on;plot(Predict(11).crosscoAC(700:1050)/max(Predict(11).crosscoAC(700:1050)),'b','LineWidth',2)
                % %         hold on;plot(Predict(11).LearnedNorm(700:1050),'k','LineWidth',2)
                % %         xlim([0 350]); ylim([0 1.1])
                % %     subplot(224);plot(Predict(7).crosscoAC(340:690)/max(Predict(7).crosscoAC(340:690)),'b','LineWidth',2)
                % %         hold on;plot(prediction11_7(700:1050)/max(prediction11_7(700:1050)),'r','LineWidth',2)
                % %         hold on;plot(Predict(7).LearnedNorm(340:690),'k','LineWidth',2)
                % %         xlim([0 350]); ylim([0 1.1])
                % %     % 11 has been cleaned; these clearly demonstrate a flatness
                % %     % early for 7 relative to 11.
                % %         subplot(221);plot(Predict(11).ResidAC(700:1050,[2:7 9:10 12:14 16:24]))
                % %         subplot(223);plot(Predict(7).ResidAC(340:690,1:20))
                  % ALSO include DC - make the MAIN POINT that the predictions are excellent for BOTH.
                  % Make the secondary, point that there is a non-significant trend
                  % for subtleties in the residual structure to mirror subtleties in
                  % learning structure.
% Here is old contingsim analysis
                % % ContingSim - BestScoreSelf2
                % % Metric 1- BestScoreSelf
                %       g=[1:22];
                %      for i=1:22
                %          ad=find(g~=i);
                % 
                %          indic2=find(BestSimX(ad,i)>0);
                %          BestScoreSelf(i)=(BestSimX(i,i));
                %          BestScoreOther(i)=mean(BestSimX(ad(indic2),i));
                %      end
                %        
                %      g=[1:22];
                %      for i=1:22
                %          ad=find(g~=i);
                % 
                %          indic2=find(Score(ad,i)>0);
                %          BestScoreSelfCS(i)=(Score(i,i));
                %          BestScoreOtherCS(i)=mean(Score(ad(indic2),i));
                %      end           
                %                 
                % % % %             
                % % % %             %figure;plot(mean(abc7_15))
                % % % %             
                % % % %         % A3Bi - 
                % % % %             figure;plot(Predict(1
                % % % %         % A3Di - experiment 1 comparison
                % % % %             % 
                % % % %             prediction7_15=jccrossco(Predict(7).ResidAC,Predict(15).Targeting-410);
                % % % %             
                % % % %             figure;
                % % % %             
                % % % % %             abc15=ContingSim2(Predict(15).Targeting-32,Predict(15).ResidAC,70);
                % % % %             hold on;plot(mean(abc15(:,750:1100))/max(mean(abc15(:,750:1100))),'b')
                % % % %             hold on;plot(Predict(15).LearnedNorm(750:1100),'k')
                % % % % %             abc15_1=ContingSim2(Predict(15).Targeting-532-32,Predict(1).ResidAC,70);
                % % % %             abc15_7=ContingSim2(Predict(15).Targeting-410-32,Predict(7).ResidAC,70);
                % % % %             hold on;plot(mean(abc15_1(:,218:568))/0.0186,'r')
                % % % % %             abc15_19=ContingSim2(Predict(15).Targeting-450-32,Predict(19).ResidAC,70);
                % % % %             hold on;plot(mean(abc15_19(:,300:650))/0.013,'g')
                % % % %             
                % % % %             figure;
                % % % %             %abc1=ContingSim2(Predict(1).Targeting-32,Predict(1).ResidAC,70);
                % % % %             hold on;plot(mean(abc1(:,218:568))/max(mean(abc1(:,218:568))),'b')
                % % % %             hold on;plot(Predict(1).LearnedNorm(218:568),'k')
                % % % %             abc1_15=ContingSim2(Predict(1).Targeting+532-32,Predict(15).ResidAC,70);
                % % % %             hold on;plot(mean(abc1_15(:,750:1100))/max(mean(abc1_15(:,750:1100))),'r')
                % % % %             %abc15_19=ContingSim2(Predict(15).Targeting-450-32,Predict(19).ResidAC,70);
                % % % %             hold on;plot(mean(abc15_19(:,300:650))/0.013,'g')
                % % % % 
                % % % %             
                % % % %             figure;plot(Predict(19).LearnedNorm(300:650),'k')
                % % % % %             abc19_15=ContingSim2(Predict(19).Targeting+450-32,Predict(15).ResidAC,70);
                % % % % %             abc19=ContingSim2(Predict(19).Targeting-32,Predict(19).ResidAC,70);
                % % % %             hold on;plot(mean(abc19_15(:,750:1100))/max(mean(abc19_15(:,750:1100))),'g')
                % % % %             hold on;plot(mean(abc19(:,300:650))/max(mean(abc19(:,300:650))),'b')
                % % % %             
                % % % %             figure;plot(Predict(1).LearnedNorm)
                % % % %             hold on;plot(Predict(1).crosscoAC/max(Predict(1).crosscoAC),'r')
                % % % %             %prediction15_1=jccrossco(Predict(15).ResidAC,Predict(1).Targeting+532);
                % % % %             hold on;plot(prediction15_1(532:end)/max(prediction15_1),'k')      
                % % % %             xlim([218 568]);ylim([0 1.1])
                % % % %         % A3Dii - experiment 15 comparison
                % % % %             figure;plot(Predict(15).LearnedNorm(750:1100));
                % % % %             hold on;plot(Predict(15).crosscoAC(750:1100)/max(Predict(15).crosscoAC(750:1100)),'r')
                % % % %             %prediction1_15=jccrossco(Predict(1).ResidAC,Predict(15).Targeting-532);
                % % % %             hold on;plot(prediction1_15(218:562)/max(prediction1_15(218:562)),'k')
                % % % %             
                % % % %             
                % % % %             for i=1:28
                % % %                 middle=Predict(i).onset+0.5*(Predict(i).offset-Predict(i).onset);
                % % %                 if isequal(Predict(i).direction,'down')
                % % %                     prcval=40;
                % % %                 else
                % % %                     prcval=60;
                % % %                 end
                % % %                 abc=ContingSim2(middle,Predict(i).ResidAC,70);
                % % %                 tt(i,:)=abc(middle-150:middle+150)./max(abc(middle-150:middle+150));
                % % %             end
                % % %             mmo=tt(:,70)+tt(:,230);
% Figure A5
            indLes=[7 13 14 15 16 18 19 21 22 25 26 27 28];
            % Make the CS predictions for INA long notes 
                for j=indLes(1:9)
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    if isequal('up',Predict(j).direction);value=70;else value=30;end
                    for i=j
                       NewTargs=Predict(i).onset+NormTargs;
                       ind1=find(NewTargs>1);
                       ind2=find(NewTargs(ind1)<1800);
                       CSLMAN(i).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidINA,value)));
                    end
                end
            % Make the CS predictions for INA short notes 
                for j=indLes(10:13)
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    if isequal('up',Predict(j).direction);value=70;else value=30;end
                    for i=j
                        NewTargs=Predict(i).onset+NormTargs;
                        ind1=find(NewTargs>1);
                        ind2=find(NewTargs(ind1)<1800);
                        CSLMAN(i).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidINA,value)));
                    end
                end
            % Fit the CS predictions for long notes - optimal scaling factor
                for i=indLes(1:9)
                    for j=i
                        count=0;
                        for k=0.5:0.01:2
                            count=count+1;
                            choices(count)=sum(abs((k*abs(CSLMAN(i).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSLMAN(i).data(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                        end
                        disterLMAN(i)=min(choices);
                    end
                end
            % Fit the CS predictions for short notes - optimal scaling factor
                for i=indLes(10:13)
                    for j=i
                        count=0;
                        for k=0.5:0.01:2
                            count=count+1;
                            choices(count)=sum(abs((k*abs(CSLMAN(i).data(Predict(i).onset:Predict(i).onset+160))./max(abs(CSLMAN(i).data(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                        end
                        disterLMAN(i)=min(choices);
                    end
                end           
            % normalize scores
                BestLMANCSnorm(1:22)=(disterLMAN(1:22)/350);
                BestLMANCSnorm(23:28)=(disterLMAN(23:28)/160);
             [h,k]=ttest(BestLMANCSnorm(indLes)-BestSelfCSnorm(indLes)) 
              % p=0.001

            %%% Get FWHM data
                    inas=indLes;
                    for ii=1:length(inas)
                        k=inas(ii);
                        middle=Predict(k).onset+(Predict(k).offset-Predict(k).onset)/2;
                        clear ac
                        clear ina
                        for i=1:size(Predict(k).ResidAC,2)
                            ac(i,:)=((xcorr(Predict(k).ResidAC(middle-80:middle+80,i))));
                        end
                        for i=1:size(Predict(k).ResidINA,2)
                            ina(i,:)=((xcorr(Predict(k).ResidINA(middle-80:middle+80,i))));
                        end
                        ccac(k).data=mean(ac)./max(mean(ac));
                        ccina(k).data=mean(ina)./max(mean(ina));
                        dat(ii)=sum(ccac(k).data-ccina(k).data)/length(ccac(k).data);
                    end
                    for ii=indLes
                        indiAC=min(find(ccac(ii).data>0.5));
                        fwhmAC(ii)=160-indiAC;
                        indiINA=min(find(ccina(ii).data>0.5));
                        fwhmINA(ii)=160-indiINA;
                    end
            % Check CVratios - why we exclude #3,#10,#24 from indLes
                % 20% cut-off
                for ii=[3 7 10 13 14 15 16 18 19 21 22 24:28]
                    CVratio(ii)=1-mean(std(Predict(ii).ResidINA(Predict(ii).onset+50:Predict(ii).offset,:)'))./mean(std(Predict(ii).ResidAC(Predict(ii).onset+50:Predict(ii).offset,:)'));
                end
            
                %                 for i=indLes
                %                     middle=Predict(i).onset+0.5*(Predict(i).offset-Predict(i).onset);
                %                     abc=(jccrossco2(Predict(i).ResidINA,middle));
                %                     hhLES(i,:)=abc(middle-150:middle+150)./max(abc(middle-150:middle+150));
                %                 end
                %                 parameter1LES=hhLES(:,70)+hhLES(:,170);

            %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%
            t=0:1:110;
            indLes=[7 13 14 15 16 18 19 21 22 25 26 27 28];
            %%% FIGURE 4 - LMAN influence %%%%%%%%%
            figure;hold on;
            % Fig.4A - example
            subplot(241)
            plot(Predict(7).ResidAC(Predict(7).onset:Predict(7).onset+350,1:20))
            xlim([0 350]);ylim([-0.06 0.06])
            subplot(242)
            plot(Predict(7).ResidINA(Predict(7).onset:Predict(7).onset+350,1:20))
            xlim([0 350]);ylim([-0.06 0.06])
            subplot(245)
            plot(Predict(27).ResidAC(Predict(27).onset:Predict(27).onset+160,41:60))
            xlim([0 160]);ylim([-0.06 0.06])
            subplot(246)
            plot(Predict(27).ResidINA(Predict(27).onset:Predict(27).onset+160,70:89))
            xlim([0 160]);ylim([-0.06 0.06])
                        %             subplot(3,3,3);hold on;
                        %             % LMAN INAs
                        %             plot(((fwhmINA(indLes([1:3 10:11]))*2)/8),((fwhmAC(indLes([1:3 10:11]))*2)/8),'*','Color','r')
                        %             % LMAN Lesions 
                        %             plot(((fwhmINA(indLes(4:9))*2)/8),((fwhmAC(indLes(4:9))*2)/8),'*','Color','b')
                        %             % RA AP-5
                        %             plot(((fwhmINA(indLes(12:13))*2)/8),((fwhmAC(indLes(12:13))*2)/8),'*','Color','g')
                        %             % How to make marker size give indication of SE???
                        %             plot(mean(fwhmINA(indLes)/4),mean(fwhmAC(indLes)/4),'+','MarkerSize',15,'Color','k')
                        %             plot(t,t,'-')
                        %             xlim([0 20]);ylim([0 20])

            subplot(243)
            % Fits %%% 7 AC - 1.03 %%% 7 LMAN - 1 %%% 27 AC - 1.05 %%% 27 LMAN - 1.09
            plot(Predict(7).LearnedNorm(Predict(7).onset:Predict(7).onset+350),'k','LineWidth',2)
            hold on;plot(abs(CSs2(7,7).data(Predict(7).onset:Predict(7).onset+350))/max(abs(CSs2(7,7).data(Predict(7).onset:Predict(7).onset+350)))*1.03,'b','LineWidth',2)
            plot(abs(CSLMAN(7).data(Predict(7).onset:Predict(7).onset+350))/max(abs(CSLMAN(7).data(Predict(7).onset:Predict(7).onset+350)))*1,'r','LineWidth',2)
            t=0:1:350;
            plot(t,mean(Predict(7).LearnedNorm(Predict(7).onset:Predict(7).onset+350)),'g','LineWidth',2)
            xlim([0 350]);ylim([0 1.1])
            subplot(247)
            plot(Predict(27).LearnedNorm(Predict(27).onset:Predict(27).onset+160),'k','LineWidth',2)
            hold on;plot(abs(CSs2(27,27).data(Predict(27).onset:Predict(27).onset+160))/max(abs(CSs2(27,27).data(Predict(27).onset:Predict(27).onset+160)))*1.05,'b','LineWidth',2)
            plot(abs(CSLMAN(27).data(Predict(27).onset:Predict(27).onset+160))/max(abs(CSLMAN(27).data(Predict(27).onset:Predict(27).onset+160)))*1.09,'r','LineWidth',2)
            t=0:1:350;
            plot(t,mean(Predict(27).LearnedNorm(Predict(27).onset:Predict(27).onset+160)),'g','LineWidth',2)
            xlim([0 160]);ylim([0 1.1])
            subplot(244)
            hold on;plot(BestSelfCSnorm(indLes([1:3 10:11])),BestLMANCSnorm(indLes([1:3 10:11])),'.','MarkerSize',20,'Color','r')
            plot(BestSelfCSnorm(indLes(4:9)),BestLMANCSnorm(indLes(4:9)),'.','MarkerSize',20,'Color','b')
            plot(BestSelfCSnorm(indLes(12:13)),BestLMANCSnorm(indLes(12:13)),'.','MarkerSize',20,'Color','g')
            hold on;plot(t,t,'-','Color','k')
            xlim([0 0.2]);ylim([0 0.2])
            plot(mean(BestSelfCSnorm(indLes)),mean(BestLMANCSnorm(indLes)),'+','MarkerSize',15,'Color','k')

            subplot(248)
            hold on;plot(4,BestSelfCSnorm(indLes),'.','MarkerSize',20,'Color','b')
                plot(3,BestOtherCSnorm(indLes),'.','MarkerSize',20,'Color','k')
                plot(3,mean(BestOtherCSnorm(indLes)),'.','MarkerSize',50,'Color','k')
                hold on;plot([3 4],[BestOtherCSnorm(indLes);BestSelfCSnorm(indLes)],'-','Color','k')
                plot(4,mean(BestSelfCSnorm(indLes)),'.','MarkerSize',50,'Color','k')
                hold on;plot([2 3],[BestLMANCSnorm(indLes);BestOtherCSnorm(indLes)],'-','Color','k')
                hold on;plot(2,[BestLMANCSnorm(indLes)],'.','MarkerSize',20,'Color','r')
                plot(2,median(BestLMANCSnorm(indLes)),'.','MarkerSize',50,'Color','k')
                hold on;plot(1,[BestScoreDCnorm(indLes)],'.','MarkerSize',20,'Color','g')
                plot(1,median(BestScoreDCnorm(indLes)),'.','MarkerSize',50,'Color','k')
                hold on;plot([1 2],[BestScoreDCnorm(indLes);BestLMANCSnorm(indLes)],'-','Color','k')
                xlim([0.7 4.3]);ylim([0 0.25])
          % Simplified 7.19
            figure;hold on;
                plot(2,BestSelfCSnorm(indLes),'.','MarkerSize',20,'Color','b')
                plot(2,mean(BestSelfCSnorm(indLes)),'.','MarkerSize',50,'Color','k')
                hold on;plot([1 2],[BestLMANCSnorm(indLes);BestSelfCSnorm(indLes)],'-','Color','k')
                hold on;plot(1,[BestLMANCSnorm(indLes)],'.','MarkerSize',20,'Color','r')
                plot(1,mean(BestLMANCSnorm(indLes)),'.','MarkerSize',50,'Color','k')
                ylim([0 0.25])
                xlim([0.5 2.5])
            

            %%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%
            % This demonstrates how DMP has a fluctuating correlational
            % structure, whereas DMP+AFP has a steadily decaying correlational
            % structure.
                for i=indLes
                start=Predict(i).onset;
                abc=(jccrossco2(Predict(i).ResidINA,start+50));
                hh2LES(i,:)=abc(start+50:start+150)./max(abc(start+50:start+150));
                end
                for i=indLes
                start=Predict(i).onset;
                abc=(jccrossco2(Predict(i).ResidAC,start+50));
                hh2(i,:)=abc(start+50:start+150)./max(abc(start+50:start+150));
                end
                figure;plot(hh2LES','r')
                hold on;plot(hh2','k')
% Figure A4
            % From experiment 7
                %         for i=1:size(bk50ht,2)
                %             bk50htprc(:,i)=(bk50e7(:,i)'-mean(bk50e7(:,1:444)'))./mean(bk50e7(:,1:444)');
                %         end
                % Change ColorMap to scaled and limits from min=0 and max=0.03
                    for i=1:size(DShifts(2).pitchALL,2)
                        dsprc2(:,i)=(DShifts(2).pitchALL(:,i)'-mean(DShifts(2).pitchALL(:,1:152)'))./mean(DShifts(2).pitchALL(:,1:152)');
                    end
                    createfigure4B(DShifts(1).pitchPRC(DShifts(1).onset:DShifts(1).offset,:)');syn
                    xlim([0 400]);ylim([0 1200])
                    createfigure4B(DShifts(2).pitchPRC(DShifts(2).onset:DShifts(2).offset,:)');syn
                    xlim([0 400]);ylim([0 size(DShifts(2).pitchALL,2)])
        %             zscores2B=(DShifts(2).pitchALL(round(median(DShifts(2).toffset))+192,:)-mean(DShifts(2).pitchALL(round(median(DShifts(2).toffset))+192,1:152)))/std(DShifts(2).pitchALL(round(median(DShifts(2).toffset))+192,1:152));
                    createfigure4B(bk50htprc(Predict(4).onset+7:Predict(4).offset,:)');syn
                    xlim([0 400]);ylim([0 size(bk50htprc,2)])

                    % Calculate predictions
                    indALL1=ContingSimAB2(DShifts(1).toffset-32,DShifts(1).pitchBaseline,DShifts(1).dirA,DShifts(1).dirB);
                    DShifts(1).CSpred=mean(DShifts(1).ResidBase(:,indALL1)');
                    indALL2=ContingSimAB2(DShifts(2).toffset-32,DShifts(2).pitchBaseline,DShifts(2).dirA,DShifts(2).dirB);
                    DShifts(2).CSpred=mean(DShifts(2).ResidBase(:,indALL2)');
                    for i=1:132
                        nums(i)=length(find(indALL1==i));
                    end
                    indReal1=find(nums>50);
                    for i=1:152
                        nums(i)=length(find(indALL2==i));
                    end
                    indReal2=find(nums>50);

                for i=1:size(pitchALL2,2)
                     or92prct2(:,i)=(pitchALL2(:,i)'-mean(pitchALL2(:,1:282)'))./mean(pitchALL2(:,1:282)');
                end
                    zscores2B=(DShifts(2).pitchALL(round(median(DShifts(2).toffset))+192,:)-mean(DShifts(2).pitchALL(round(median(DShifts(2).toffset))+192,1:152)))/std(DShifts(2).pitchALL(round(median(DShifts(2).toffset))+192,1:152));
                    figure;hold on


        %%%% CALCULATIONS
        % DShifts data set
                    trackjA=zeros(6,2000);
                    trackk1=[];
                    for i=1:6
                        j=((DShifts(i).pitchALL(-32+round(median(DShifts(i).toffset)),size(DShifts(i).pitchBaseline,2):end)-mean(DShifts(i).pitchBaseline(-32+round(median(DShifts(i).toffset)),:))));
                        trackjA(i,1:length(j))=j;
                        if isequal(DShifts(i).dirA,'up')
                            fac=1;
                        else
                            fac=-1;
                        end

                        for ii=1:20
                            trackk1(i,ii)=median(trackjA(i,ii*70-69:ii*70))*fac;
                        end
                    end
                    trackjB=zeros(6,2000);
                    trackk2=[];
                    for i=1:6
                        j=((DShifts(i).pitchALL(-32+round(median(DShifts(i).toffset))+192,size(DShifts(i).pitchBaseline,2):end)-mean(DShifts(i).pitchBaseline(-32+round(median(DShifts(i).toffset))+192,:))));
                        trackjB(i,1:length(j))=j;
                        if isequal(DShifts(i).dirB,'up')
                            fac=1;
                        else
                            fac=-1;
                        end
                        for ii=1:20
                            trackk2(i,ii)=median(trackjB(i,ii*70-69:ii*70))*fac;
                        end
                    end

        %             % y-axis=magnitude of adaptive shift
        %                 figure;plot(trackk2')
        %                 hold on;plot(trackk2','.','MarkerSize',20)
        %                 plot(trackk1')
        %                 plot(trackk1','.','MarkerSize',20,'Color','k')

           %%% SUMMARY PLOT             
                        for i=1:6
                            if isequal(DShifts(i).dirA,'up')
                                facA=1;
                            else
                                facA=-1;
                            end
                            numlong=size(DShifts(i).pitchALL,2)-size(DShifts(i).pitchBaseline,2);
                            nummid=round(numlong/2);
                            endptA(i)=-1*facA*median(trackjA(i,numlong-100:numlong));
                            endptB(i)=-1*facA*median(trackjB(i,numlong-100:numlong));
                            midptA(i)=-1*facA*median(trackjA(i,nummid-50:nummid+50));
                            midptB(i)=-1*facA*median(trackjB(i,nummid-50:nummid+50));
                            startptA(i)=0;
                            startptB(i)=0;
                        end
                        figure;hold on;
                         plot(1,startptB,'*')
                         plot(2,midptB,'.','MarkerSize',20)
                         plot(3,endptB,'.','MarkerSize',20)
                         plot(1,startptA,'*')
                         plot(2,midptA,'.','MarkerSize',20)
                         plot(3,endptA,'.','MarkerSize',20)
                         plot([1 2 3],[startptB;midptB;endptB],'b')
                         plot([1 2 3],[startptA;midptA;endptA],'r')





        %%%% THE MAIN FIGURE EXAMPLE
                    figure;subplot(3,2,3);hold on;
                    plot(DShifts(1).ResidBase,'k')
                    plot(DShifts(1).ResidBase(:,indReal1),'r')
                    plot(DShifts(1).CSpred,'b','LineWidth',2)
                    plot([DShifts(1).onset DShifts(1).offset],[0 0],'g')
                    [xS,yS]=hist(DShifts(1).toffset-32,20);
                    stairs(yS,-0.05+0.2*(xS/length(DShifts(1).toffset)),'LineWidth',2,'Color','b')
                    stairs(yS+192,-0.05+0.2*(xS/length(DShifts(1).toffset)),'LineWidth',2,'Color','b')
                    xlim([DShifts(1).onset DShifts(1).offset])
                    ylim([-0.05 0.05])
                    subplot(3,2,5);hold on;
                    plot(DShifts(2).ResidBase,'k')
                    plot(DShifts(2).ResidBase(:,indReal2),'r')
                    plot(DShifts(2).CSpred,'b','LineWidth',2)
                    plot([DShifts(2).onset DShifts(2).offset],[0 0],'g')
                    [xS,yS]=hist(DShifts(2).toffset-32,20);
                    stairs(yS,-0.05+0.1*(xS/length(DShifts(2).toffset)),'LineWidth',2,'Color','b')
                    stairs(yS+192,-0.05+0.1*(xS/length(DShifts(2).toffset)),'LineWidth',2,'Color','b')
                    xlim([DShifts(2).onset DShifts(2).offset])
                    ylim([-0.05 0.05])
                    subplot(3,2,1);hold on;
                    mtar=mean(Predict(4).Targeting)
                    mprc=prctile(Predict(4).ResidAC(round(mtar),1:150),70)
                    indH=find(Predict(4).ResidAC(round(mtar),1:150)>mprc)
                    plot(Predict(4).ResidAC,'k')
                    hold on;plot(Predict(4).ResidAC(:,indH),'r')
                    hold on;plot(mean(Predict(4).ResidAC(:,indH)'),'b','LineWidth',2)   
                    plot([Predict(4).onset Predict(4).offset],[0 0],'g')
                    [xS,yS]=hist(Predict(4).Targeting-32,20);
                    stairs(yS,-0.05+0.2*(xS/length(Predict(4).Targeting)),'LineWidth',2,'Color','b')
                    xlim([Predict(4).onset+7 Predict(4).offset])
                    ylim([-0.05 0.05])

                    subplot(3,2,2);hold on;
                    plot(Predict(4).LearnedNorm,'k')
                    plot(abs(CSs64(4,4).data)/max(abs(CSs64(4,4).data(Predict(4).onset+7:Predict(4).offset)))*1.05)
                    [xS,yS]=hist(Predict(4).Targeting-32,20);
                    stairs(yS,0+2*(xS/length(Predict(4).Targeting)),'LineWidth',2,'Color','b')
                    xlim([Predict(4).onset+7 Predict(4).offset]);ylim([0 1.1])
                    subplot(3,2,4);hold on;
                    plot((mean(DShifts(1).pitchALL(:,700:end)')-mean(DShifts(1).pitchBaseline'))/max(abs(mean(DShifts(1).pitchALL(DShifts(1).onset:DShifts(1).offset,700:end)')-mean(DShifts(1).pitchBaseline(DShifts(1).onset:DShifts(1).offset,:)'))),'k')
                    plot(DShifts(1).CSpred*80)
                    plot([DShifts(1).onset DShifts(1).offset],[0 0],'g')
                    [xS,yS]=hist(DShifts(1).toffset-32,20);
                    stairs(yS,-1.3+2*(xS/length(DShifts(1).toffset)),'LineWidth',2,'Color','b')
                    stairs(yS+192,-1.3+2*(xS/length(DShifts(1).toffset)),'LineWidth',2,'Color','b')
                    xlim([DShifts(1).onset DShifts(1).offset]);ylim([-1.3 1])
                    subplot(3,2,6);hold on;
                    plot((mean(DShifts(2).pitchALL(:,500:end)')-mean(DShifts(2).pitchBaseline'))/max(abs(mean(DShifts(2).pitchALL(DShifts(2).onset:DShifts(2).offset,500:end)')-mean(DShifts(2).pitchBaseline(DShifts(2).onset:DShifts(2).offset,:)'))),'k')
                    plot(DShifts(2).CSpred*80)
                    plot([DShifts(2).onset DShifts(2).offset],[0 0],'g')
                    [xS,yS]=hist(DShifts(2).toffset-32,20);
                    stairs(yS,-1.3+1*(xS/length(DShifts(2).toffset)),'LineWidth',2,'Color','b')
                    stairs(yS+192,-1.3+1*(xS/length(DShifts(2).toffset)),'LineWidth',2,'Color','b')
                    xlim([DShifts(2).onset DShifts(2).offset]);ylim([-1.3 1])
        %%% CV curve example
            for i=1:28
                if ~isempty(Predict(i).ResidINA)
                    CVs(i)=1-mean(std(Predict(i).ResidINA(Predict(i).onset:Predict(i).offset,:)'))./mean(std(Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:)'));
                end
            end
            figure;plot(1,CVs(indLes([1:3 10:11])),'+','MarkerSize',10,'Color','b')
            hold on;plot(2,CVs(indLes([4:9])),'+','MarkerSize',10,'Color','r')
            plot(3,CVs(indLes([12:13])),'+','MarkerSize',10,'Color','g')
            plot(1,CVs(indLes([1:3 10:11])),'.','MarkerSize',15,'Color','b')
            hold on;plot(2,CVs(indLes([4:9])),'.','MarkerSize',15,'Color','r')
            plot(3,CVs(indLes([12:13])),'.','MarkerSize',15,'Color','g')
            xlim([0.7 3.3])
            ylim([0 0.6])



        %%% A4H - structure of adaptation - first blue line is first
        %%% contingency point, second blue line is second contingency point
            figure;hold on
            for i=1:length(DShifts)
                if isequal(DShifts(i).dirB,'up')
                    factor=1;
                else
                    factor=-1;
                end
                shifted(i).data=(mean(DShifts(i).pitchALL(:,end-200:end)')-mean(DShifts(i).pitchBaseline'))
                ons=-192+(DShifts(i).onset-median(DShifts(i).toffset));
                x=ons:1:ons+(DShifts(i).offset-DShifts(i).onset);
                if factor==1
                    plot(x,factor*shifted(i).data(DShifts(i).onset:DShifts(i).offset),'Color','k','LineWidth',2)
                end
                if factor==-1
                    plot(x,factor*shifted(i).data(DShifts(i).onset:DShifts(i).offset),'Color','r','LineWidth',2)
                end
            end
            plot([0-32 0-32],[-40 50],'b')
            plot([-192-32 -192-32],[-40 50],'b')
            plot([-400 300],[0 0],'b')
            xlim([-370 280])
            
       %%%% A4H2 - Michael's figure with the maximal adaptive changes.
       figure;hold on
            for i=1:length(DShifts)
                ons=-96+(DShifts(i).onset-median(DShifts(i).toffset));
                [a,b]=max(shifted(i).data(DShifts(i).onset:DShifts(i).offset));
                [c,d]=min(shifted(i).data(DShifts(i).onset:DShifts(i).offset));
                if isequal(DShifts(i).dirB,'up')
                    plot(b+ons,a,'.','MarkerSize',20,'Color','k'); plot(d+ons,c,'.','MarkerSize',20,'Color','k')
                    plot([b+ons b+ons],[0 a],'k');plot([d+ons d+ons],[0 c],'k')
                    
                else
                    i
                    plot(b+ons,a,'.','MarkerSize',20,'Color','r'); plot(d+ons,c,'.','MarkerSize',20,'Color','r')
                    plot([b+ons b+ons],[0 a],'r');plot([d+ons d+ons],[0 c],'r')
                end
            end
            plot([0-32 0-32],[-100 100]);plot([-300 250],[0 0])
            plot([92-32 92-32],[-100 100],'g');plot([-92-32 -92-32],[-100 100],'g')
            ylim([-80 80]);xlim([-300 200])
        %%% A4G - how good a fit
                for i=1:8
                    CSpred(i).data=ContingSimAB2(DShifts(i).toffset-32,DShifts(i).pitchBaseline,DShifts(i).dirA,DShifts(i).dirB,DShifts(i).onset,DShifts(i).offset);

                    if isequal(DShifts(i).dirB,'up')
                        prc=51;
                    else
                        prc=49;
                    end
                    CSpredB(i).data=mean(ContingSim2B(DShifts(i).toffset+192-32,jc_residuals(DShifts(i).pitchBaseline),prc,DShifts(i).dirB));

                end
                
                
            %%%% SENSITIVITY
                for i=8
                    CSpred46(i).data=ContingSimAB3(DShifts(i).toffset-32,DShifts(i).pitchBaseline,DShifts(i).dirA,DShifts(i).dirB,DShifts(i).onset,DShifts(i).offset);
                end

                
            

%                  figure;
%                 for i=1:6
%                      hold on;subplot(2,3,i);hold on;
%                     if isequal(DShifts(i).dirB,'up')
%                         factor=1;
%                     else
%                         factor=-1;
%                     end
%                     facB=abs(shifted(i).data(round(median(DShifts(i).toffset))+192-32))/abs(CSpredB(i).data(round(median(DShifts(i).toffset))+192-32));
%                     fac=abs(shifted(i).data(round(median(DShifts(i).toffset))+192-32))/abs(CSpred(i).data(round(median(DShifts(i).toffset))+192-32));
% 
%                     prAB=factor*CSpred(i).data(DShifts(i).onset:DShifts(i).offset)*fac;
%                     prB=factor*CSpredB(i).data(DShifts(i).onset:DShifts(i).offset)*facB;
%                     Actual=factor*shifted(i).data(DShifts(i).onset:DShifts(i).offset);
%                      plot(prAB/max(abs(Actual)),'r')
%                      plot(Actual/max(abs(Actual)))
%                      plot(prB/max(abs(Actual)),'k')
%                     distanceAB(i)=sum(abs(prAB-Actual))/(DShifts(i).offset-DShifts(i).onset+1);
%                     distanceAB(i)=distanceAB(i)/max(abs(Actual));
%                     distanceB(i)=sum(abs(prB-Actual))/(DShifts(i).offset-DShifts(i).onset+1);
%                     distanceB(i)=distanceB(i)/max(abs(Actual));
%                 end
                % Optimally normalize
                figure;
                for i=1:8
                     hold on;subplot(2,3,i);hold on;
                    if isequal(DShifts(i).dirB,'up')
                        factor=1;
                    else
                        factor=-1;
                    end
                    facB=max(abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset)))/max(abs(CSpredB(i).data(DShifts(i).onset:DShifts(i).offset)));
                    facAB=max(abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset)))/max(abs(CSpred(i).data(DShifts(i).onset:DShifts(i).offset)));

                    Actual=factor*shifted(i).data(DShifts(i).onset:DShifts(i).offset);
                    TTab=[0:facAB/50:facAB*2];
                    for j=1:length(TTb)
                        prAB(j).data=factor*CSpred(i).data(DShifts(i).onset:DShifts(i).offset)*TTab(j);
                        distanceABx(j)=sum(abs(prAB(j).data-Actual))/(DShifts(i).offset-DShifts(i).onset+1);
                    end
                    [aa,bb1]=min(distanceABx);
                    plot(prAB(bb1).data/max(abs(Actual)),'r')
                    TTb=[0:facB/50:facB*2];
                    for j=1:length(TTb)
                        prB(j).data=factor*CSpredB(i).data(DShifts(i).onset:DShifts(i).offset)*TTb(j);
                        distanceBx(j)=sum(abs(prB(j).data-Actual))/(DShifts(i).offset-DShifts(i).onset+1);
                    end
                    [aa,bb2]=min(distanceBx);
                    plot(prB(bb2).data/max(abs(Actual)),'k')
                    
                     plot(Actual/max(abs(Actual)))

                     distanceAB(i)=sum(abs(prAB(bb1).data-Actual))/(DShifts(i).offset-DShifts(i).onset+1);
                     distanceAB(i)=distanceAB(i)/max(abs(Actual));
                     distanceB(i)=sum(abs(prB(bb2).data-Actual))/(DShifts(i).offset-DShifts(i).onset+1);
                     distanceB(i)=distanceB(i)/max(abs(Actual));
                end

                    

                     actual=zeros(8,2000);
                    predicted=zeros(8,2000);

                    for i=1:length(DShifts)
                        if isequal(DShifts(i).dirB,'up')
                            factor=1;
                        else
                            factor=-1;
                        end
                        shifted(i).data=(mean(DShifts(i).pitchALL(:,end-200:end)')-mean(DShifts(i).pitchBaseline'))
                        ons=-192+(DShifts(i).onset-median(DShifts(i).toffset));
                        x=ons:1:ons+(DShifts(i).offset-DShifts(i).onset);
                        maximum=max(abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset)));
                        newx=x+400;
                        actual(i,round(newx))=factor*(1/maximum)*shifted(i).data(DShifts(i).onset:DShifts(i).offset);
                        predicted(i,round(newx))=factor*CSpred(i).data(DShifts(i).onset:DShifts(i).offset);
                        predicted(i,round(newx))=predicted(i,round(newx))/max(abs((predicted(i,round(newx)))));
                        plot(x,predicted(i,round(newx)),'r')
                        plot(x,actual(i,round(newx)),'Color','k','LineWidth',2)
                    end
                    %hold on;plot([-400:1:1599],bjf,'b')
                    plot([0-32 0-32],[-1.1 1.1],'b')
                    plot([-192-32 -192-32],[-1.1 1.1],'b')
                    plot([-400 300],[0 0],'b')
                    xlim([-330 120])
                    ylim([-1.1 1.1])
                    for i=1:2000
                        mnactual(i)=0;
                        mnpredicted(i)=0;
                        ind=find(actual(:,i)~=0)
                        mnactual(i)=mean(actual(ind,i));
                        mnpredicted(i)=mean(predicted(ind,i));
                    end
                    hold on;plot([-400:1:1599],mnactual,'LineWidth',4,'Color','k')
                    hold on;plot([-400:1:1599],mnpredicted,'LineWidth',4,'Color','r')
            



                figure;hold on;
                plot([1 2],[distanceB;distanceAB],'Color','k')
                plot(1,distanceB,'.','MarkerSize',20,'Color','r')
                plot(2,distanceAB,'.','MarkerSize',20,'Color','b')
                xlim([0.8 2.2]);ylim([0 0.7])
                
            %%%%
            % 7.16.09
            %%%
  
      %%% This figure shows the raw data for predictions and actual and indicates that 
      %%% the mean predictions and mean actuals look very similar.
      
      
      %%%%%% FIGURE 4F
                    figure;hold on
                    actual=zeros(8,2000);
                    predicted=zeros(8,2000);

                    for i=1:length(DShifts)
                        if isequal(DShifts(i).dirB,'up')
                            factor=1;
                        else
                            factor=-1;
                        end
                        shifted(i).data=(mean(DShifts(i).pitchALL(:,end-200:end)')-mean(DShifts(i).pitchBaseline'))
                        ons=-192+(DShifts(i).onset-median(DShifts(i).toffset));
                        x=ons:1:ons+(DShifts(i).offset-DShifts(i).onset);
                        maximum=max(abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset)));
                        newx=x+400;
                        actual(i,round(newx))=factor*(1/maximum)*shifted(i).data(DShifts(i).onset:DShifts(i).offset);
                        predicted(i,round(newx))=factor*CSpred46(i).data(DShifts(i).onset:DShifts(i).offset);
                        predicted(i,round(newx))=predicted(i,round(newx))/max(abs((predicted(i,round(newx)))));
                        plot(x,predicted(i,round(newx)),'r')
                        plot(x,actual(i,round(newx)),'Color','k','LineWidth',2)
                    end
                    %hold on;plot([-400:1:1599],bjf,'b')
                    plot([0-32 0-32],[-1.1 1.1],'b')
                    plot([-192-32 -192-32],[-1.1 1.1],'b')
                    plot([-400 300],[0 0],'b')
                    xlim([-330 120])
                    ylim([-1.1 1.1])
                    for i=1:2000
                        mnactual(i)=0;
                        mnpredicted(i)=0;
                        ind=find(actual(:,i)~=0)
                        mnactual(i)=mean(actual(ind,i));
                        mnpredicted(i)=mean(predicted(ind,i));
                    end
                    hold on;plot([-400:1:1599],mnactual,'LineWidth',4,'Color','k')
                    hold on;plot([-400:1:1599],mnpredicted,'LineWidth',4,'Color','r')
            
            
            %%% DEcay in 20ms
                                    % Actual data from single contingency
                                    % experiments
                                              figure;hold on
                                              for i=1:28
                                                  [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                                                  abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28
                                                     b=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                                                     b=round(b);
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+b);
                                                     abb(i,700:700+dister2)=Predict(i).LearnedNorm(Predict(i).onset+b:Predict(i).offset);
                                                 end
                                                 count=0;
                                                 for i=1:28
                                                     if abb(i,540)~=0 && abb(i,860)~=0
                                                         count=count+1;
                                                         onecont20(count)=mean([abb(i,540) abb(i,860)]);
                                                     else if abb(i,860)~=0
                                                         count=count+1;
                                                         onecont20(count)=abb(i,860);
                                                         else if abb(i,540)~=0
                                                                 count=count+1;
                                                                 onecont20(count)=abb(i,540);
                                                             end
                                                         end
                                                     end
                                                 end
                                                 
                                                 
                                % Predictions from single contingency
                                % experiments
                                              for i=1:28
                                                  aax=CSs64(i,i).data;
                                                  if isequal(Predict(i).direction,'up')
                                                      a(i)=max(aax(Predict(i).onset:Predict(i).offset));
                                                  else
                                                      a(i)=min(aax(Predict(i).onset:Predict(i).offset));
                                                  end

                                                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                                                  abb=aax(Predict(i).onset:Predict(i).onset+btop);
                                                  left=length(abb);
                                                  abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                                                  right=length(abb)-left;
                                                  t=-1*left:1:right-1;
                                                  abb=abb/a(i);
                                                  plot(t/8,abb,'Linewidth',2,'Color','k')
                                              end
                                                 for i=1:28
                                                     notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                                                 end
                                                 abb=zeros(28,1400);
                                                 for i=1:28
                                                     b=median(Predict(i).Targeting)-Predict(i).onset; %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                                                     dister1=b;
                                                     dister2=notewidth(i)*8-b;

                                                     abb(i,700-dister1:700)=abs(CSs64(i,i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                                                     abb(i,700:700+dister2)=abs(CSs64(i,i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                                                 end
                                                 count=0;
                                                 for i=1:28
                                                     if abb(i,540)~=0 && abb(i,860)~=0
                                                         count=count+1;
                                                         onecont20pred(count)=mean([abb(i,540) abb(i,860)]);
                                                     else if abb(i,860)~=0
                                                         count=count+1;
                                                         onecont20pred(count)=abb(i,860);
                                                         else if abb(i,540)~=0
                                                                 count=count+1;
                                                                 onecont20pred(count)=abb(i,540);
                                                             end
                                                         end
                                                     end
                                                 end

                                                 
                                                 
                    %% FIGURE 4H                             
                          figure;plot(onecont20pred,onecont20,'.','MarkerSize',20,'Color','r')  
                          hold on;plot(mean(onecont20pred),mean(onecont20),'+','MarkerSize',15,'Color','k')
                          plot(predicted(:,240-32),actual(:,240-32),'.','MarkerSize',20,'Color','k')    
                          plot(mean(predicted(:,240-32)),mean(actual(:,240-32)),'+','MarkerSize',15,'Color','k')
                          xlim([-1 1])
                          ylim([-1 1])
                          
                          

             % Learning is driven by the contingencies: shape --- p=0.00076 compared with zero shift
                 for i=1:5
                     distDC(i)=mean((abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset))/max(abs(shifted(i).data(DShifts(i).onset:DShifts(i).offset)))));
                 end
                 [h,p]=ttest(dist,distDC)
                          
                          
                          
% FIGURE A6
        figure;hold on
        for i=1:9%length(indLes)
            ind=indLes(i);
            middle=round((Predict(ind).offset-Predict(ind).onset)/2+Predict(ind).onset);
            [x,psv(i).data]=jcpsd2(Predict(ind).ResidAC(middle-150:middle+150,:),8000);
            [x,psvINA(i).data]=jcpsd2(Predict(ind).ResidINA(middle-150:middle+150,:),8000);
            plot(x,median(psv(i).data))
            plot(x,median(psvINA(i).data),'r')
        end
        for i=1:9
            bff=median(psv(i).data)./median(psvINA(i).data);
            bfrats(i,:)=bff(1:10);
        end
        for i=1:6
            middle=round((AlldataZFlesion(i).offset-AlldataZFlesion(i).onset)/2+AlldataZFlesion(i).onset);
            [x,ZFpsv(i).data]=jcpsd2(AlldataZFlesion(i).pitchUDpre(middle-150:middle+150,:),8000);
            [x,ZFpsvINA(i).data]=jcpsd2(AlldataZFlesion(i).pitchUDpost(middle-150:middle+150,:),8000);
        end
        for i=1:6
            zff=median(ZFpsv(i).data)./median(ZFpsvINA(i).data);
            zfrats(i,:)=zff(1:10);
        end
        for i=1:9
            bff=median(psv(i).data)./median(psvINA(i).data);
            bfrats(i,:)=bff(1:10);
        end
        figure;hold on;subplot(121);plot(x(1:10),bfrats)
        xlim([10 100])
        subplot(122);plot(x(1:10),zfrats)
        xlim([10 100])
%%%
            % Make the CS predictions for long notes 
                for j=1:22
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    if isequal('up',Predict(j).direction);value=51;else value=49;end
                    for i=1:22
                       NewTargs=Predict(i).onset+NormTargs;
                       ind1=find(NewTargs>1);
                       ind2=find(NewTargs(ind1)<1800);
                       CSs2(i,j).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidAC,value)));
                    end
                end
            % Make the CS predictions for short notes 
                for j=23:28
                    j
                    NormTargs=Predict(j).Targeting-Predict(j).onset; % dist b/f offset
                    if isequal('up',Predict(j).direction);value=51;else value=49;end
                    for i=23:28
                        NewTargs=Predict(i).onset+NormTargs;
                        ind1=find(NewTargs>1);
                        ind2=find(NewTargs(ind1)<1800);
                        CSs2(i,j).data=(mean(ContingSim2(NewTargs(ind1(ind2)),Predict(i).ResidAC,value)));
                    end
                end
            % Fit the CS predictions for long notes - optimal scaling factor
                for i=1:22
                    for j=1:22
                        count=0;
                        for k=0.5:0.01:2
                            count=count+1;
                            choices(count)=sum(abs((k*abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+350))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+350)));
                        end
                        disters2(i,j)=min(choices);
                    end
                end
            % Fit the CS predictions for short notes - optimal scaling factor
                for i=23:28
                    for j=23:28
                        count=0;
                        for k=0.5:0.01:2
                            count=count+1;
                            choices(count)=sum(abs((k*abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+160))./max(abs(CSs2(i,j).data(Predict(i).onset:Predict(i).onset+160))))-Predict(j).LearnedNorm(Predict(j).onset:Predict(j).onset+160)));
                        end
                        disters2(i,j)=min(choices);
                    end
                end           
            % Plot the *unscaled* CS predictions against the actual learning
                for i=1:28
                    subplot(5,6,i);hold on;plot(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).onset+350))./max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).onset+350))))
                    plot(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350),'r')
                    xlim([0 350]);ylim([0 1.1])
                end
            % DC - how good is this prediction
                    % DC - goodness of fit 
                    for i=1:22
                        BestScoreDC(i)=sum(abs(mean(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350))-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+350)));
                    end
                    for i=23:28
                        BestScoreDC(i)=sum(abs(mean(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+160))-Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+160)));
                    end                 
            % Self/Other - how good are these predictions
                for i=1:28
                    BestSelfCS(i)=disters2(i,i);
                end
                % long - other
                g=[1:22];
                for i=1:22
                    indn=find(g~=i);
                    BestOtherCS(i)=mean(disters2(indn,i));
                end
                % short - other
                g=[23:28];
                for i=23:28
                    indn=find(g~=i);
                    BestOtherCS(i)=mean(disters2(g(indn),i));
                end
                
                % normalize scores
                    BestSelfCSnorm(1:22)=(BestSelfCS(1:22)/350);
                    BestSelfCSnorm(23:28)=(BestSelfCS(23:28)/160);
                    BestOtherCSnorm(1:22)=(BestOtherCS(1:22)/350);
                    BestOtherCSnorm(23:28)=(BestOtherCS(23:28)/160);
                    BestScoreDCnorm(1:22)=(BestScoreDC(1:22)/350);
                    BestScoreDCnorm(23:28)=(BestScoreDC(23:28)/160);
%%%%
%%%%%%
%%%%%%

% beast.m can calculate as a proportion, but this actually reduces the
% difference, so not worthwhile

% average sigma of adaptation
for i=1:28
    sigmashift(i)=max(Predict(i).Learned(Predict(i).onset:Predict(i).offset))/(median(std(Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:)'))*mmean(i));
end
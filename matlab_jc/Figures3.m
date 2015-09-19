% Figure 1

% The experiment here is with bk50w18


% Figure 1A
    % Spectrogram for or92 End_Up experiment
    % The way to do it is to label an early note 'b' and then used tafsimjc
    % to take the 'average' spectrogram (but it's just that one note).
    % Use or92 instead of bk50 because in bk50 the white noise is notched.
    
    % Before: OR92FFbaseline=[mean(pitch1117C(680:740,:)),mean(pitch1118A(680:740,:)),mean(pitch1118B(680:740,:))];
    % FF distribution before (see Figure 1B)
        [x1,y1]=hist(OR92FFbaseline);
        x1(1)=0;
        x1(end)=0;
        figure;hold on;
        stairs(x1/max(x1),y1,'b','LineWidth',2)
        plot([0 1],[2350 2350],'k','LineWidth',2)
% Figure 1B
    % FF distributions before and after
    % Before: OR92FFbaseline=[mean(pitch1117C(680:740,:)),mean(pitch1118A(680:740,:)),mean(pitch1118B(680:740,:))];
    % After: OR92FFshifted=[mean(pitch1123A(680:740,:)),mean(pitch1123B(680:740,:)),mean(pitch1123C(680:740,:)),mean(pitch1124A(680:740,:))];
        [x1,y1]=hist(OR92FFbaseline);
        [x2,y2]=hist(OR92FFshifted);
        x1(1)=0;
        x1(end)=0;
        x2(1)=0;
        x2(end)=0;
        figure;hold on;
        stairs(x1/max(x1),y1,'b','LineWidth',2)
        stairs(x2/max(x2),y2,'r','LineWidth',2)
% Figure 2
    % Sample targeting histogram is the only thing I need here
    % PredTarg6=Predict(6).Targeting
    [x,y]=hist(PredTarg6(find(PredTarg6<700))); 
    x(end)=0;x(1)=0;
    figure;stairs(y,x,'k','LineWidth',2)
    hold on;plot(510,400,'.','MarkerSize',20,'Color','k')
    ylim([0 500])
% Figure 3A - spectrograms pre and post
    % Spectrograms of or92or82 upshift before and after with 1st harmonic
    % FF contour.  
      % AD7=Alldatafirstten1sig(7);
        % createfigureA1A.m
        t=BK50Base.t;
        f=BK50Base.f;
        distance=112; % 30ms-16ms - distance between onset of spec and onset of contour
        initpt=370;
        finpt=870;
        meanBase=(mean(AD7.baselineAC'));
        meanShifted=(0.5*(mean(AD7.exp(17).selectedpitchcurves')+mean(AD7.exp(18).selectedpitchcurves')));
        
        figure;subplot(121);imagesc(t,f,log(BK50Base.av));syn;ylim([0 10000]);xlim([0.04 0.14])
        hold on;plot(t1(initpt-distance:finpt-distance)+0.032,meanBase(initpt:finpt),'LineWidth',2)
        subplot(122);imagesc(t,f,log(BK50Shift.av));syn;ylim([0 10000]);xlim([0.04 0.14])
        hold on;plot(t1(initpt-distance:finpt-distance)+0.032,meanShifted(initpt:finpt),'LineWidth',2)
        [x,y]=hist(Predict(4).Targeting-Predict(4).onset,20);
        x=[0 x];
        y=[y(1) y];
        stairs(t1(round(y+initpt))+0.032,x*50,'LineWidth',2)
        % adjust to min=5 and max=15 and colormap=hot
        
% Figure 3B - one syllable, 3 shifts

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
          plot([median(Predict(4).Targeting)-Predict(4).onset median(Predict(4).Targeting)-Predict(4).onset],[ 0.3],'.','MarkerSize',15,'Color','b')
          plot([median(Predict(5).Targeting)-Predict(5).onset median(Predict(5).Targeting)-Predict(5).onset],[ 0.3],'.','MarkerSize',15,'Color','k')
          plot([median(Predict(6).Targeting)-Predict(6).onset median(Predict(6).Targeting)-Predict(6).onset],[ 0.3],'.','MarkerSize',15,'Color','r')
          [a,b4]=max(Predict(4).LearnedNorm(Predict(4).onset:Predict(4).offset));
          [a,b5]=max(Predict(5).LearnedNorm(Predict(5).onset:Predict(5).offset));
          [a,b6]=max(Predict(6).LearnedNorm(Predict(6).onset:Predict(6).offset));
          plot([b4 b4],[ 1.1],'.','MarkerSize',15,'Color','b')
          plot([b5 b5],[ 1.1],'.','MarkerSize',15,'Color','k')
          plot([b6 b6],[ 1.1],'.','MarkerSize',15,'Color','r')
% Figure 3C - All learning curves centered around median targeting position

            figure;hold all
            for i=[1:28]
                [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                left=length(abb);
                abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                right=length(abb)-left;
                t=-1*left:1:right-1;
                plot(t/8,abb,'Linewidth',2)
%                 if i==30 || i==31
%                     plot(t/8,abb,'Linewidth',2,'Color','g')
%                 end
            end
            for i=1:28
                notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
            end
            abb=zeros(28,1400);
            for i=[1:28]
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
                ind=find(abb([1:28],i)>0);
                if ~isempty(ind)
                    mnabb(i)=mean(abb(ind,i));
                    seabb(i)=std(abb(ind,i))/sqrt(length(ind));
                end
            end
            t=-542:1:559;
            hold on;plot(t/8,mnabb(158:1259),'r','Linewidth',3)
            ylim([0 1.05])
            xlim([-25 25])


% Figure 3D - median targeting position vs. maximal FF change

    % This code allows stdev bars but doesn't really add anything.
            %         figure;hold on;
            %         for i=1:28
            %            [a,tops(i)]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
            %            middle(i)=(median(Predict(i).Targeting)-Predict(i).onset)/8;
            %            middleplus(i)=(median(Predict(i).Targeting)+std(Predict(i).Targeting)-Predict(i).onset)/8;
            %            middleminus(i)=(median(Predict(i).Targeting)-std(Predict(i).Targeting)-Predict(i).onset)/8;
            %            tops(i)=tops(i)/8;
            %            plot([middleminus(i) middle(i) middleplus(i)],[tops(i) tops(i) tops(i)],'-','Color','k')
            %         end
        for i=1:28
           [a,tops(i)]=max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
           middle(i)=median(Predict(i).Targeting)-Predict(i).onset;
        end
            % get in milliseconds 
        middle=middle/8;
        tops=tops/8;
        figure;hold on;
        plot(middle,tops,'.','MarkerSize',15,'Color','k')
        plot(middleplus,tops,'.','MarkerSize',15)
        p=polyfit(middle,tops,1);
        t=0:0.1:100;
        hold on;plot(t,p(2)+p(1)*t,'b')
        xlim([0 70])
        ylim([0 70])
% Figure 3E - Model 1 is wrong
        counter=zeros(1,6);
        for i=1:28
            edgeleft=round(median(Predict(i).Targeting)-2*std(Predict(i).Targeting));
            edgeright=round(median(Predict(i).Targeting)+2*std(Predict(i).Targeting));
            if edgeleft>Predict(i).onset
                counter(1)=counter(1)+1;
                neg0(counter(1))=Predict(i).LearnedNorm(edgeleft);
                if edgeleft-80>Predict(i).onset
                    counter(2)=counter(2)+1;
                    neg10(counter(2))=Predict(i).LearnedNorm(edgeleft-80);
                    if edgeleft-160>Predict(i).onset
                        counter(3)=counter(3)+1;
                        neg20(counter(3))=Predict(i).LearnedNorm(edgeleft-160);
                    end
                end
            end
            if edgeright<Predict(i).offset
                counter(4)=counter(4)+1;
                pos0(counter(4))=Predict(i).LearnedNorm(edgeright);
                if edgeright+80<Predict(i).offset
                    counter(5)=counter(5)+1;
                    pos10(counter(5))=Predict(i).LearnedNorm(edgeright+80);
                    if edgeright+160<Predict(i).offset
                        counter(6)=counter(6)+1;
                        pos20(counter(6))=Predict(i).LearnedNorm(edgeright+160);
                    end
                end
            end
        end
        counter=zeros(1,6);
        for i=30:31
            edgeleft=round(median(Predict(i).Targeting)-2*std(Predict(i).Targeting));
            edgeright=round(median(Predict(i).Targeting)+2*std(Predict(i).Targeting));
            if edgeleft>Predict(i).onset
                counter(1)=counter(1)+1;
                ZFneg0(counter(1))=Predict(i).LearnedNorm(edgeleft);
                if edgeleft-80>Predict(i).onset
                    counter(2)=counter(2)+1;
                    ZFneg10(counter(2))=Predict(i).LearnedNorm(edgeleft-80);
                    if edgeleft-160>Predict(i).onset
                        counter(3)=counter(3)+1;
                        ZFneg20(counter(3))=Predict(i).LearnedNorm(edgeleft-160);
                    end
                end
            end
            if edgeright<Predict(i).offset
                counter(4)=counter(4)+1;
                ZFpos0(counter(4))=Predict(i).LearnedNorm(edgeright);
                if edgeright+80<Predict(i).offset
                    counter(5)=counter(5)+1;
                    ZFpos10(counter(5))=Predict(i).LearnedNorm(edgeright+80);
                    if edgeright+160<Predict(i).offset
                        counter(6)=counter(6)+1;
                        ZFpos20(counter(6))=Predict(i).LearnedNorm(edgeright+160);
                    end
                end
            end
        end        
        figure;plot(0,1,'.','MarkerSize',20,'Color','b')
        hold on;plot(-4.5,neg20,'.','MarkerSize',20,'Color','b')
        plot(-3,neg10,'.','MarkerSize',20,'Color','b')
        plot(-2,neg0,'.','MarkerSize',20,'Color','b')
        plot(2,pos0,'.','MarkerSize',20,'Color','b')
        plot(3,pos10,'.','MarkerSize',20,'Color','b')
        plot(4.5,pos20,'.','MarkerSize',20,'Color','b')
        % incorporate ZF experiments
        plot(-4.5,ZFneg20,'.','MarkerSize',20,'Color','r')
        plot(-3,ZFneg10,'.','MarkerSize',20,'Color','r')
        plot(-2,ZFneg0,'.','MarkerSize',20,'Color','r')
        plot(2,ZFpos0,'.','MarkerSize',20,'Color','r')
        plot(3,ZFpos10,'.','MarkerSize',20,'Color','r')
        plot(4.5,ZFpos20,'.','MarkerSize',20,'Color','r')
        plot([-4.5 -3 -2 0 2 3 4.5],[mean([neg20 ZFneg20]) mean([neg10 ZFneg10]) mean([neg0 ZFneg0]) 1 mean([pos0 ZFpos0]) mean([pos10 ZFpos10]) mean([pos20 ZFpos20])],'k')
        g=gaussian([1:600],300,200);
        t=[-3:1/100:3];
        plot(t(1:end-1),g)
        
        
% Figure 3F - shows how predictions are generated
            figure;hold on;
            %             subplot(241);plot(Predict(7).ResidAC(Predict(7).onset+30:Predict(7).onset+350,70:90))
            %             xlim([0 320]);ylim([-0.06 0.06])
            mtar=mean(Predict(7).Targeting)
            mprc=prctile(Predict(7).ResidAC(round(mtar),:),60)
            indH=find(Predict(7).ResidAC(round(mtar),41:end)>mprc)
            subplot(242);hold on;plot(Predict(7).ResidAC(:,41:end),'k')

            hold on;plot(Predict(7).ResidAC(:,40+indH),'b')
            hold on;plot(CSs64(7,7).data,'r','LineWidth',4)
            xlim([Predict(7).onset+30 Predict(7).onset+350])
            ylim([-0.06 0.05])
            [xS,yS]=hist(Predict(7).Targeting-32,20);
            stairs(yS+30,-0.06+0.2*(xS/(2*length(Predict(7).Targeting))),'LineWidth',4,'Color','g')
            plot([0 1000],[0 0],'y','LineWidth',2)

            subplot(243);plot(Predict(7).LearnedNorm(Predict(7).onset+30:Predict(7).onset+350),'k','LineWidth',2)
            hold on;plot(abs(CSs64(7,7).data(Predict(7).onset+30:Predict(7).onset+350))/max(abs(CSs64(7,7).data(Predict(7).onset+30:Predict(7).onset+350)))*1.00,'r','LineWidth',2)
            stairs(yS-Predict(7).onset,xS/sum(xS),'LineWidth',4,'Color','g')
            t=0:1:350;
            xlim([0 320]);ylim([0 1.1])

        
% Figure 3G - compare predicted vs. actual -
% RUN CODE for 3C first!!!!!!!
                % CS predictions - variability with targeting - centered at targ position
                for i=1:28
                    aax=CS73(i).data;

                        a(i)=max(aax(Predict(i).onset:Predict(i).offset));

                    btop=median(Predict(i).Targeting)-Predict(i).onset;
                    abb=aax(Predict(i).onset:Predict(i).onset+btop);
                    left=length(abb);
                    abb=[abb aax(Predict(i).onset+btop:Predict(i).offset)];
                    right=length(abb)-left;
                    t=-1*left:1:right-1;
                    abb=abb/a(i);
                end
                for i=1:28
                    notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                end
                abb=zeros(28,1400);
                for i=1:28
                    b=round(median(Predict(i).Targeting)-Predict(i).onset); %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                    dister1=(b);
                    dister2=(notewidth(i)*8-b);
                    abb(i,700-dister1:700)=abs(CS73(i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                    abb(i,700:700+dister2)=abs(CS73(i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
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
                %%%%%
                
                
                
                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                figure;hold on;
                plot([t/8;t/8],[mnabb(158:1259)/j1+seabb(158:1259)/j1;mnabb(158:1259)/j1-seabb(158:1259)/j1],'color','k')
                plot(t/8,mnabbT(158:1259)/j2,'r','LineWidth',3) % targ imprecision included
                xlim([-40 40]);ylim([0 1.05])
                plot(t/8,mnabb(158:1259)/j1)
                xlim([-25 25])
              % How much of the shape is explained by the predictions?
                % 25ms on either side
                    % r=0.9694
                    % r^2=0.9397
                
                
% Figure 4 - inter-syllable comparison

% #7 and #24 - all this bullshit about #30 is ridiculous....
figure;hold on;
plot(CSs64(7,7).data(373:553)/0.0175)
plot(CSs64(24,24).data(243:423)/0.017/0.99/0.922)
plot(Predict(7).LearnedNorm(373:553),'r')
plot(Predict(24).LearnedNorm(243:423),'r')
ylim([0 1.05])
xlim([0 180])






% Figure 4 - lines ~2300 in 416.m
            figure;subplot(221);plot(Predict(1).LearnedNorm(Predict(1).onset+30:Predict(1).onset+350),'k','LineWidth',2)
                hold on;plot(abs(CSs64(1,1).data(Predict(1).onset+30:Predict(1).onset+350))/max(abs(CSs64(1,1).data(Predict(1).onset+30:Predict(1).onset+350)))*1.00,'r','LineWidth',2)
                plot(abs(CSs64(4,1).data(Predict(4).onset+30:Predict(4).onset+350))/max(abs(CSs64(4,1).data(Predict(4).onset+30:Predict(4).onset+350)))*1.01,'b','LineWidth',2)   
                t=0:1:350;
                plot(t,mean(Predict(1).LearnedNorm(Predict(1).onset+30:Predict(1).onset+350)),'g','LineWidth',2)
                xlim([0 320]);ylim([0 1.1])
                [xS,yS]=hist(Predict(1).Targeting-32,20);
                stairs(yS-Predict(1).onset,xS/sum(xS),'LineWidth',4,'Color','k')

            subplot(223);plot(Predict(4).LearnedNorm(Predict(4).onset+30:Predict(4).onset+350),'k','LineWidth',2)
                hold on;plot(abs(CSs64(4,4).data(Predict(4).onset+30:Predict(4).onset+350))/max(abs(CSs64(4,4).data(Predict(4).onset+30:Predict(4).onset+350)))*1.05,'r','LineWidth',2)
                plot(abs(CSs64(1,4).data(Predict(1).onset+30:Predict(1).onset+350))/max(abs(CSs64(1,4).data(Predict(1).onset+30:Predict(1).onset+350)))*0.99,'b','LineWidth',2)   
                t=0:1:350;
                plot(t,mean(Predict(4).LearnedNorm(Predict(4).onset:Predict(4).onset+350)),'g','LineWidth',2)
                xlim([0 320]);ylim([0 1.1]) 
                [xS,yS]=hist(Predict(4).Targeting-32,20);
                stairs(yS-Predict(4).onset,xS/sum(xS),'LineWidth',4,'Color','k')

            subplot(222);hold on;
            plot(2,BestSelfCSnorm64,'.','MarkerSize',20,'Color','r')
                plot(2,mean(BestSelfCSnorm64),'.','MarkerSize',50,'Color','k')
                hold on;plot([1 2],[BestOtherCSnorm64;BestSelfCSnorm64],'-','Color','k')
                hold on;plot(1,[BestOtherCSnorm64],'.','MarkerSize',20,'Color','b')
                plot(1,mean(BestOtherCSnorm64),'.','MarkerSize',50,'Color','k')
                hold on;plot(3,[BestScoreDCnorm64],'.','MarkerSize',20,'Color','g')
                plot(3,mean(BestScoreDCnorm64),'.','MarkerSize',50,'Color','k')
                hold on;plot([2 3],[BestSelfCSnorm64;BestScoreDCnorm64],'-','Color','k')
                xlim([0.7 3.3]);ylim([0 0.25])
                
% Figure "5" - adjacent syllables                
       % cardinal/Adjsylls901.mat
       
       % Figure 5A
          figure;hold on;subplot(211);
          imagesc(Exp12.t,Exp12.f,log(Exp12.av));ylim([0 10000]);xlim([-0.01 0.14]);syn
          subplot(212);hold on;plot([750:800],pred12(750:800),'r','Linewidth',3)
          plot([750:800],(mean(Adjsylls(12).pitchApost(750:800,:)')-mean(Adjsylls(12).pitchApre(750:800,:)'))/139.4,'k','Linewidth',3)
          plot([1150:1600],pred12(1150:1600),'r','Linewidth',3)
          plot([1150:1600],(mean(Adjsylls(12).pitchApost(1150:1600,:)')-mean(Adjsylls(12).pitchApre(1150:1600,:)'))/139.4,'k','Linewidth',3)
          xlim([713-80 1913-80])
          hold on;plot([750:800],std(Adjsylls(12).pitchApre(750:800,:)')/139.4,'Linewidth',2)
          hold on;plot([1150:1600],std(Adjsylls(12).pitchApre(1150:1600,:)')/139.4,'Linewidth',2)
          %hold on;plot([750:800],-1*std(Adjsylls(12).pitchApre(750:800,:)')/139.4,'Linewidth',2)
          %hold on;plot([1150:1600],-1*std(Adjsylls(12).pitchApre(1150:1600,:)')/139.4,'Linewidth',2)
          hold on;plot([713-80 1913-80],[0 0],'k')
          [x,y]=hist(Predict(16).Targeting);
          x=[0 x 0];
          y=[y(1) y y(end)];
          stairs(y+540,x/sum(x),'LineWidth',2,'Color','g')
          ylim([-0.5 1.1])
      % Figure 5B
      figure;plot(zscorepred,zscore,'.','Markersize',20)
      hold on;plot([0 0],[-1 1],'k')
      hold on;plot([-1 2],[0 0],'k')
      xlim([-0.5 2]);ylim([-0.5 0.5])
      
%                       %
%                           for i=[1:8 10:15]
%                               en=Adjsylls(i).expnum;
%                               if isequal(Predict(en),'up')
%                                   prct=60;
%                               else
%                                   prct=40;
%                               end
%                               CSI=ContingSimIND(mean(Adjsylls(i).regA),jc_residuals(Adjsylls(i).pitchApre),prct);
%                               sA=mean(mean(Adjsylls(i).pitchApre(Adjsylls(i).regA(1):Adjsylls(i).regA(2),CSI)'))-mean(mean(Adjsylls(i).pitchApre(Adjsylls(i).regA(1):Adjsylls(i).regA(2),:)'));
%                               sB=mean(mean(Adjsylls(i).pitchBpre(Adjsylls(i).regB(1):Adjsylls(i).regB(2),CSI)'))-mean(mean(Adjsylls(i).pitchBpre(Adjsylls(i).regB(1):Adjsylls(i).regB(2),:)'));
%                               pred(i)=sB/sA;
%                               actual(i)=Adjsylls(i).Bshift/Adjsylls(i).Ashift;
%                           end
% 
%                           for i=[1:8 10:15]
%                               if Adjsylls(i).Ashift>0
%                                   coef=1;
%                               else
%                                   coef=-1;
%                               end
%                               zscore(i)=(coef*Adjsylls(i).Bshift)/mean(std(Adjsylls(i).pitchApre(Adjsylls(i).regB(1):Adjsylls(i).regB(2),:)'));
%                           end
% 
%                           % In each of the cases where the pred's are significant, the
%                           % pitch correlations are also significant.
%                           for i=[1:8 10:15]
%                               if Adjsylls(i).Ashift>0
%                                   coef=1;
%                               else
%                                   coef=0;
%                               end
%                               zscorepred(i)=(coef*pred(i)*Adjsylls(i).Ashift)/mean(std(Adjsylls(i).pitchApre(Adjsylls(i).regB(1):Adjsylls(i).regB(2),:)'));
%                           end
          
% Figure 6 - Double contingency experiments
            %% Predictions
                for i=1:8
                    CSpred(i).data=ContingSimAB2(DShifts(i).toffset-32,jc_residuals(DShifts(i).pitchBaseline),DShifts(i).dirA,DShifts(i).dirB,DShifts(i).onset,DShifts(i).offset);
                    if isequal(DShifts(i).dirB,'up')
                        prc=40;
                    else
                        prc=60;
                    end
                end
                
                
            %% Actual
            for i=1:length(DShifts)
                if isequal(DShifts(i).dirB,'up')
                    factor=1;
                else
                    factor=-1;
                end
                shifted(i).data=(mean(DShifts(i).pitchALL(:,end-200:end)')-mean(DShifts(i).pitchBaseline'))
                ons=-192+(DShifts(i).onset-median(DShifts(i).toffset));
                x=ons:1:ons+(DShifts(i).offset-DShifts(i).onset);
            end
% Figure 6A - spectrogram showing A and B
            % bk50w18 - same as in 3A
            pointA=median(Predict(4).Targeting)/12.5+16;
            pointB=pointA+24;
            figure;subplot(121);imagesc(t,f,BK50Base.av);syn;ylim([0 10000]);xlim([0 0.15])
            hold on;plot([pointA/1000 pointA/1000],[2000 3000],'g','Linewidth',2)
            plot([pointB/1000 pointB/1000],[2000 3000],'g','Linewidth',2)
            % colormap hot min=-500 max=10000
            subplot(122);imagesc(t,f,BK50Shift.av);syn;ylim([0 10000]);xlim([0 0.2])
% Figure 6B - pred/actual
    figure;hold on;
    subplot(131);hold on;
    targpos=median(Predict(12).Targeting)-32;
    plot(CSs64(12,12).data/max(abs(CSs64(12,12).data(Predict(12).onset:Predict(12).offset))),'r','Linewidth',2)
    plot(Predict(12).LearnedNorm,'k','Linewidth',2)
    ylim([-1.1 1.1])
    plot([Predict(12).onset+7 Predict(12).offset+80],[0 0])
    plot([targpos targpos],[0.5 0.8],'LineWidth',2,'Color','g')
    xlim([targpos-240 targpos+100])

    subplot(132);hold on;
        targpos=median(DShifts(3).toffset)-32+192;
    plot(CSpred(3).data/max(abs(CSpred(3).data(DShifts(3).onset:DShifts(3).offset))),'r','Linewidth',2)
    plot(shifted(3).data/max(abs(shifted(3).data(DShifts(3).onset:DShifts(3).offset))),'k','Linewidth',2)
    ylim([-1.1 1.1])
    %[xS,yS]=hist(DShifts(3).toffset-32,20);
    plot([targpos-192 targpos-192],[0.5 0.8],'LineWidth',2,'Color','g')
    plot([targpos targpos],[0.5 0.8],'LineWidth',2,'Color','g')
    plot([DShifts(3).onset DShifts(3).offset],[0 0])
    xlim([targpos-240 targpos+100])
 %%%%%%%
    % Determine the amount of decay within 12ms of median targeting
    % Decay=[(FF at WN onset)-(FF 12ms before WN onset)]/[maximal adaptive change in FF]
    % Choose 12ms - halfway between
    
    % NEEDS MORE WORK - WHAT IS THE APPROPRIATE METRIC
  % Single contingency experiments
    width=96; % 12ms - allows n=28
    prepoint=[700-width];
    whichind=abb(:,prepoint)~=0;
    decaySingle=[(abb(whichind,700))-(abb(whichind,prepoint))]/1;

  % Double contingency experiments
  for i=1:8
      WNpos=round(median(DShifts(i).toffset)-32+192);
      coef=1-2*isequal(DShifts(i).dirB,'down');
      absshift=coef*shifted(i).data;
      normshifts(i).data=absshift/max(absshift(DShifts(i).onset:DShifts(i).offset));
      decayDouble(i)=[normshifts(i).data(WNpos)-normshifts(i).data(WNpos-width)]/1;
  end
figure;plot([0.9:0.01:1.05],1-decaySingle,'V')
hold on;plot([0.93:0.01:1.00],1-decayDouble,'V','Color','r')

%%%%%%
% Center all curves on the time of WN delivery
% Normalize all curves to maximal adaptive shift in the B direction
        t=[100:-1:-300];
        clear predicted
        clear actual
            for i=1:8
                WNpos=round(median(DShifts(i).toffset)-32+192);
                coef=1-2*isequal(DShifts(i).dirB,'down');
                absshift=coef*shifted(i).data;
                normshifts(i).data=absshift/max(absshift(DShifts(i).onset:DShifts(i).offset));
                actual(i,:)=normshifts(i).data(WNpos+t);
                %%
                abspredshift=coef*CSpred(i).data;
                prednormshifts(i).data=abspredshift/max(abspredshift(DShifts(i).onset:DShifts(i).offset));
                predicted(i,:)=prednormshifts(i).data(WNpos+t);
            end
%% load /cardinal/FiguresA/Alldata216.mat
            figure;hold on;
           
            scaling=mean(actual(:,101));
            scalingP=mean(predicted(:,101));
            plot([t;t],[mean(actual)/scaling-(std(actual)/sqrt(8))/scaling;mean(actual)/scaling+(std(actual)/sqrt(8))/scaling],'k')  
            plot([t;t],[mean(predicted)/scalingP-(std(predicted)/sqrt(8))/scalingP;mean(predicted)/scalingP+(std(predicted)/sqrt(8))/scalingP],'r')
%             ciplot(mean(actual)/scaling-(std(actual)/sqrt(8))/scaling,mean(actual)/scaling+(std(actual)/sqrt(8))/scaling,t,'k',0.8)
%             ciplot(mean(predicted)/scalingP-(std(predicted)/sqrt(8))/scalingP,mean(predicted)/scalingP+(std(predicted)/sqrt(8))/scalingP,t,'r',0.8)
            % Run code for 3C and 3G
            t2=[-400:1:559];
            plot([t2;t2],[mnabb(300:1259)/mnabb(700)-seabb(300:1259)/mnabb(700);mnabb(300:1259)/mnabb(700)+seabb(300:1259)/mnabb(700)],'g')
            plot([0 0],[-2 1.5],'k')
            plot([-192 -192],[-2 1.5],'k')
            plot([-400 300],[0 0],'k')
            xlim([-300 80])
            figure;hold on;
            plot(t,mean(actual)/scaling,'k','Linewidth',3)
            plot(t2,mnabb(300:1259)/mnabb(700),'g','Linewidth',3)
            plot([0 0],[-2 1.5],'k')
            plot([-192 -192],[-2 1.5],'k')
            plot([-400 300],[0 0],'k')
            width=100;
            decayDouble=actual(:,100+width);%./actual(:,101);
            prepoint=[700-width];
            whichind=abb(:,prepoint)~=0;
            decaySingle=abb(whichind,prepoint);%./abb(whichind,700);
            plot(-width,decayDouble,'V','Markersize',10,'Linewidth',1,'Color','k')
            plot(-width,decaySingle,'V','Markersize',10,'Linewidth',1,'Color','g')
            ylim([-2 1.5])
            xlim([-300 80])
            width=80;
            decayDouble=actual(:,100+width);%./actual(:,101);
            prepoint=[700-width];
            whichind=abb(:,prepoint)~=0;
            decaySingle=abb(whichind,prepoint);%./abb(whichind,700);
            plot(-width,decayDouble,'V','Markersize',10,'Linewidth',1,'Color','k')
            plot(-width,decaySingle,'V','Markersize',10,'Linewidth',1,'Color','g')

% for i=1:8
% plot([-300:1:100]/8,normshifts(i).data(round(median(DShifts(i).toffset)-32+192)-300:round(median(DShifts(i).toffset)-32+192)+100),'r','Linewidth',2)
% end
xlim([-20 5])
% Figure 6C - mean +/-SE
                figure;hold on;
                     actual=zeros(8,2000);% Figure 6 - Double contingency experiments
            %% Predictions
                for i=1:8
                    CSpred(i).data=ContingSimAB2(DShifts(i).toffset-32,jc_residuals(DShifts(i).pitchBaseline),DShifts(i).dirA,DShifts(i).dirB,DShifts(i).onset,DShifts(i).offset);
                    if isequal(DShifts(i).dirB,'up')
                        prc=40;
                    else
                        prc=60;
                    end
                end
            %% Actual
            for i=1:length(DShifts)
                if isequal(DShifts(i).dirB,'up')
                    factor=1;
                else
                    factor=-1;
                end
                shifted(i).data=(mean(DShifts(i).pitchALL(:,end-200:end)')-mean(DShifts(i).pitchBaseline'))
                ons=-192+(DShifts(i).onset-median(DShifts(i).toffset));
                x=ons:1:ons+(DShifts(i).offset-DShifts(i).onset);
            end
% Figure 6A - spectrogram showing A and B
            % bk50w18 - same as in 3A
            pointA=median(Predict(4).Targeting)/12.5+16;
            pointB=pointA+24;
            figure;subplot(121);imagesc(t,f,BK50Base.av);syn;ylim([0 10000]);xlim([0 0.15])
            hold on;plot([pointA/1000 pointA/1000],[2000 3000],'g','Linewidth',2)
            plot([pointB/1000 pointB/1000],[2000 3000],'g','Linewidth',2)
            % colormap hot min=-500 max=10000
            subplot(122);imagesc(t,f,BK50Shift.av);syn;ylim([0 10000]);xlim([0 0.2])
% Figure 6B - pred/actual
    figure;hold on;
    subplot(131);hold on;
    targpos=median(Predict(12).Targeting)-32;
    plot(CSs64(12,12).data/max(abs(CSs64(12,12).data(Predict(12).onset:Predict(12).offset))),'r','Linewidth',2)
    plot(Predict(12).LearnedNorm,'k','Linewidth',2)
    ylim([-1.1 1.1])
    plot([Predict(12).onset+7 Predict(12).offset+80],[0 0])
    plot([targpos targpos],[0.5 0.8],'LineWidth',2,'Color','g')
    xlim([targpos-240 targpos+100])

    subplot(132);hold on;
        targpos=median(DShifts(3).toffset)-32+192;
    plot(CSpred(3).data/max(abs(CSpred(3).data(DShifts(3).onset:DShifts(3).offset))),'r','Linewidth',2)
    plot(shifted(3).data/max(abs(shifted(3).data(DShifts(3).onset:DShifts(3).offset))),'k','Linewidth',2)
    ylim([-1.1 1.1])
    %[xS,yS]=hist(DShifts(3).toffset-32,20);
    plot([targpos-192 targpos-192],[0.5 0.8],'LineWidth',2,'Color','g')
    plot([targpos targpos],[0.5 0.8],'LineWidth',2,'Color','g')
    plot([DShifts(3).onset DShifts(3).offset],[0 0])
    xlim([targpos-240 targpos+100])
 %%%%%%%
    % Determine the amount of decay within 12ms of median targeting
    % Decay=[(FF at WN onset)-(FF 12ms before WN onset)]/[maximal adaptive change in FF]
    % Choose 12ms - halfway between
    
    % NEEDS MORE WORK - WHAT IS THE APPROPRIATE METRIC
  % Single contingency experiments
    width=96; % 12ms - allows n=28
    prepoint=[700-width];
    whichind=abb(:,prepoint)~=0;
    decaySingle=[(abb(whichind,700))-(abb(whichind,prepoint))]/1;

  % Double contingency experiments
  for i=1:8
      WNpos=round(median(DShifts(i).toffset)-32+192);
      coef=1-2*isequal(DShifts(i).dirB,'down');
      absshift=coef*shifted(i).data;
      normshifts(i).data=absshift/max(absshift(DShifts(i).onset:DShifts(i).offset));
      decayDouble(i)=[normshifts(i).data(WNpos)-normshifts(i).data(WNpos-width)]/1;
  end
figure;plot([0.9:0.01:1.05],1-decaySingle,'V')
hold on;plot([0.93:0.01:1.00],1-decayDouble,'V','Color','r')

%%%%%%
% Center all curves on the time of WN delivery
% Normalize all curves to maximal adaptive shift in the B direction
        t=[100:-1:-300];
        clear predicted
        clear actual
            for i=1:8
                WNpos=round(median(DShifts(i).toffset)-32+192);
                coef=1-2*isequal(DShifts(i).dirB,'down');
                absshift=coef*shifted(i).data;
                normshifts(i).data=absshift/max(absshift(DShifts(i).onset:DShifts(i).offset));
                actual(i,:)=normshifts(i).data(WNpos+t);
                %%
                abspredshift=coef*CSpred(i).data;
                prednormshifts(i).data=abspredshift/max(abspredshift(DShifts(i).onset:DShifts(i).offset));
                predicted(i,:)=prednormshifts(i).data(WNpos+t);
            end
%% load /cardinal/FiguresA/Alldata216.mat
            figure;hold on;
           
            scaling=mean(actual(:,101));
            scalingP=mean(predicted(:,101));
            plot([t;t],[mean(actual)/scaling-(std(actual)/sqrt(8))/scaling;mean(actual)/scaling+(std(actual)/sqrt(8))/scaling],'k')  
            plot([t;t],[mean(predicted)/scalingP-(std(predicted)/sqrt(8))/scalingP;mean(predicted)/scalingP+(std(predicted)/sqrt(8))/scalingP],'r')
%             ciplot(mean(actual)/scaling-(std(actual)/sqrt(8))/scaling,mean(actual)/scaling+(std(actual)/sqrt(8))/scaling,t,'k',0.8)
%             ciplot(mean(predicted)/scalingP-(std(predicted)/sqrt(8))/scalingP,mean(predicted)/scalingP+(std(predicted)/sqrt(8))/scalingP,t,'r',0.8)
            % Run code for 3C and 3G
            t2=[-400:1:559];
            plot([t2;t2],[mnabb(300:1259)/mnabb(700)-seabb(300:1259)/mnabb(700);mnabb(300:1259)/mnabb(700)+seabb(300:1259)/mnabb(700)],'g')
            plot([0 0],[-2 1.5],'k')
            plot([-192 -192],[-2 1.5],'k')
            plot([-400 300],[0 0],'k')
            xlim([-300 80])
            figure;hold on;
            plot(t,mean(actual)/scaling,'k','Linewidth',3)
            plot(t2,mnabb(300:1259)/mnabb(700),'g','Linewidth',3)
            plot([0 0],[-2 1.5],'k')
            plot([-192 -192],[-2 1.5],'k')
            plot([-400 300],[0 0],'k')
            width=100;
            decayDouble=actual(:,100+width);%./actual(:,101);
            prepoint=[700-width];
            whichind=abb(:,prepoint)~=0;
            decaySingle=abb(whichind,prepoint);%./abb(whichind,700);
            plot(-width,decayDouble,'V','Markersize',10,'Linewidth',1,'Color','k')
            plot(-width,decaySingle,'V','Markersize',10,'Linewidth',1,'Color','g')
            ylim([-2 1.5])
            xlim([-300 80])
            width=80;
            decayDouble=actual(:,100+width);%./actual(:,101);
            prepoint=[700-width];
            whichind=abb(:,prepoint)~=0;
            decaySingle=abb(whichind,prepoint);%./abb(whichind,700);
            plot(-width,decayDouble,'V','Markersize',10,'Linewidth',1,'Color','k')
            plot(-width,decaySingle,'V','Markersize',10,'Linewidth',1,'Color','g')

% for i=1:8
% plot([-300:1:100]/8,normshifts(i).data(round(median(DShifts(i).toffset)-32+192)-300:round(median(DShifts(i).toffset)-32+192)+100),'r','Linewidth',2)
% end
xlim([-20 5])
% Figure 6C - mean +/-SE
                figure;hold on;
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
                    end
                    %hold on;plot([-400:1:1599],bjf,'b')
                    figure;hold on;
                    plot([0-32 0-32],[-1.1 1.1],'k')
                    plot([-192-32 -192-32],[-1.1 1.1],'k')
                    plot([-400 300],[0 0],'k')
                    xlim([-330 120])
                    ylim([-1.1 1.1])
                    for i=1:2000
                        mnactual(i)=0;
                        mnpredicted(i)=0;
                        ind=find(actual(:,i)~=0)
                        mnactual(i)=mean(actual(ind,i));
                        mnSE1actual(i)=mean(actual(ind,i))+std(actual(ind,i))/sqrt(length(ind));
                        mnSE2actual(i)=mean(actual(ind,i))-std(actual(ind,i))/sqrt(length(ind));
                        mnpredicted(i)=mean(predicted(ind,i));
                        mnSE1predicted(i)=mean(predicted(ind,i))+std(predicted(ind,i))/sqrt(length(ind));
                        mnSE2predicted(i)=mean(predicted(ind,i))-std(predicted(ind,i))/sqrt(length(ind));

                    end
                    t1=[-400:1:1599];
                    plot([t1;t1],[mnSE2actual;mnSE1actual],'k')  
                    plot([t1;t1],[mnSE2predicted;mnSE1predicted],'r')
        % Run code for 3C and 3G
                t2=[-400:1:559];
                plot([t2;t2],[mnabb(300:1259)/j1-seabb(300:1259)/j1;mnabb(300:1259)/j1+seabb(300:1259)/j1],'b')

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
                    end
                    %hold on;plot([-400:1:1599],bjf,'b')
                    figure;hold on;
                    plot([0-32 0-32],[-1.1 1.1],'k')
                    plot([-192-32 -192-32],[-1.1 1.1],'k')
                    plot([-400 300],[0 0],'k')
                    xlim([-330 120])
                    ylim([-1.1 1.1])
                    for i=1:2000
                        mnactual(i)=0;
                        mnpredicted(i)=0;
                        ind=find(actual(:,i)~=0)
                        mnactual(i)=mean(actual(ind,i));
                        mnSE1actual(i)=mean(actual(ind,i))+std(actual(ind,i))/sqrt(length(ind));
                        mnSE2actual(i)=mean(actual(ind,i))-std(actual(ind,i))/sqrt(length(ind));
                        mnpredicted(i)=mean(predicted(ind,i));
                        mnSE1predicted(i)=mean(predicted(ind,i))+std(predicted(ind,i))/sqrt(length(ind));
                        mnSE2predicted(i)=mean(predicted(ind,i))-std(predicted(ind,i))/sqrt(length(ind));

                    end
                    t1=[-400:1:1599];
                    plot([t1;t1],[mnSE2actual;mnSE1actual],'k')  
                    plot([t1;t1],[mnSE2predicted;mnSE1predicted],'r')
        % Run code for 3C and 3G
                t2=[-400:1:559];
                plot([t2;t2],[mnabb(300:1259)/j1-seabb(300:1259)/j1;mnabb(300:1259)/j1+seabb(300:1259)/j1],'b')
% Figure 6D - variation - MSB figure
       figure;hold on
                   plot([92-32 92-32],[-100 100],'g');plot([-92-32 -92-32],[-100 100],'g')
            for i=1:length(DShifts)
                ons=-96+(DShifts(i).onset-median(DShifts(i).toffset));
                [a,b]=max(shifted(i).data(DShifts(i).onset:DShifts(i).offset));
                [c,d]=min(shifted(i).data(DShifts(i).onset:DShifts(i).offset));
                if isequal(DShifts(i).dirB,'up')
                    plot(b+ons,a,'.','MarkerSize',20,'Color','b'); plot(d+ons,c,'.','MarkerSize',20,'Color','b')
                    plot([b+ons b+ons],[0 a],'b');plot([d+ons d+ons],[0 c],'b')
                    plot(60,shifted(i).data(round(median(DShifts(i).toffset)+160)),'.','MarkerSize',20,'Color','b')
                    plot(-124,shifted(i).data(round(median(DShifts(i).toffset)-32)),'.','MarkerSize',20,'Color','b')
                else
                    i
                    plot(b+ons,a,'.','MarkerSize',20,'Color','r'); plot(d+ons,c,'.','MarkerSize',20,'Color','r')
                    plot([b+ons b+ons],[0 a],'r');plot([d+ons d+ons],[0 c],'r')
                    plot(60,shifted(i).data(round(median(DShifts(i).toffset)+160)),'.','MarkerSize',20,'Color','r')
                    plot(-124,shifted(i).data(round(median(DShifts(i).toffset)-32)),'.','MarkerSize',20,'Color','r')

                end
            end
            plot([0-32 0-32],[-100 100],'k');plot([-300 250],[0 0],'k')
            ylim([-80 80]);xlim([-300 200],'k')
% Figure 7 - LMAN 
% load LMANtest824.mat
% more m files in correctforn.m

%%%% Actual adaptation for indLes birds
       figure;hold on
            for i=[1:28]
                [btop]=median(Predict(i).Targeting)-Predict(i).onset;%max(Predict(i).LearnedNorm(Predict(i).onset:Predict(i).offset));
                abb=Predict(i).LearnedNorm(Predict(i).onset:Predict(i).onset+btop);
                left=length(abb);
                abb=[abb Predict(i).LearnedNorm(Predict(i).onset+btop:Predict(i).offset)];
                right=length(abb)-left;
                t=-1*left:1:right-1;
            end
            for i=[1:28]
                notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
            end
            abb=zeros(28,1400);
            for i=[1:28]
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
                ind=find(abb(indLes,i)>0);
                if ~isempty(ind)
                    mnabb(i)=mean(abb(indLes(ind),i));
                    seabb(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                end
            end
            t=-542:1:559;
            figure;hold on;
            plot([t/8;t/8;t/8],[mnabb(158:1259)-seabb(158:1259);mnabb(158:1259);mnabb(158:1259)+seabb(158:1259)],'g')
            plot(t/8,mnabb(158:1259),'k','Linewidth',3)
            ylim([0 1.05])
            xlim([-40 40])
%%% Predicted from all data for indLes birds
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
                end
                for i=1:28
                    notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                end
                abb=zeros(28,1400);
                for i=1:28
                    b=round(median(Predict(i).Targeting)-Predict(i).onset); %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                    dister1=(b);
                    dister2=(notewidth(i)*8-b);
                    abb(i,700-dister1:700)=abs(CSs64(i,i).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                    abb(i,700:700+dister2)=abs(CSs64(i,i).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                end
                mnabbT=zeros(1,1400);
                seabbT=zeros(1,1400);
                for i=1:1400
                    ind=find(abb(indLes,i)>0);
                    if ~isempty(ind)
                        mnabbT(i)=mean(abb(indLes(ind),i));
                        seabbT(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                    end
                end
                t=-542:1:559;
                %%%%%
                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                plot(t/8,mnabbT(158:1259),'r','LineWidth',3) % targ imprecision included
%%% Prediction for just DMP
        % THIS CALCULATES DATA FOR EACH EXPERIMENT
                for k=1:13
                    i=indLes(k);
                    aax=CSLMAN(k).data;
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
                end
                for i=1:28
                    notewidth(i)=(Predict(i).offset-Predict(i).onset)./8;
                end
                abb=zeros(28,1400);
                for k=1:13
                    i=indLes(k);
                    b=round(median(Predict(i).Targeting)-Predict(i).onset); %max(abs(CSs2(i,i).data(Predict(i).onset:Predict(i).offset)));
                    dister1=(b);
                    dister2=(notewidth(i)*8-b);
                    abb(i,700-dister1:700)=abs(CSLMAN(k).data(Predict(i).onset:Predict(i).onset+b)/a(i));
                    abb(i,700:700+dister2)=abs(CSLMAN(k).data(Predict(i).onset+b:Predict(i).offset)/a(i));
                end
                mnabbT=zeros(1,1400);
                seabbT=zeros(1,1400);
                for i=1:1400
                    ind=find(abb(indLes,i)>0);
                    if ~isempty(ind)
                        mnabbT(i)=mean(abb(indLes(ind),i));
                        seabbT(i)=std(abb(indLes(ind),i))/sqrt(length(indLes(ind)));
                    end
                end
                t=-542:1:559;
                %%%%%
                %%% FINAL PLOT - generates "A1H"
                j1=max(mnabb(158:1259));
                j2=max(mnabbT(158:1000));
                %j3=max(mnabbNT(158:1259));
                plot(t/8,mnabbT(158:1259),'b','LineWidth',3) % targ imprecision included

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURES FOR LMAN only component
        for i=1:13
            CVratio(indLes(i))=1-mean(std(Predictx(indLes(i)).ResidINA(Predictx(indLes(i)).onset+20:Predictx(indLes(i)).offset-20,:)'))/mean(std(Predictx(indLes(i)).ResidAC(Predictx(indLes(i)).onset+20:Predictx(indLes(i)).offset-20,:)'));
        end
        for i=1:size(mnccAC,2)
            weightedFF(:,i)=(mnccAC(:,i).^2-CVratio(indLes(i))*mnccINA(:,i).^2)/(1-CVratio(i));
            weightedFF(:,i)=weightedFF(:,i)/max(weightedFF(:,i));
        end
        for i=1:size(mnccAC,1)
            ind1=find(weightedFF(i,:)>0);
            weightedFF2(i)=mean(weightedFF(i,ind1));
        end
        weightedFF2=weightedFF2/max(weightedFF2);
        % 31.875ms (precisely)
        
        
        
        
% Figure 8D
t=0:1;
            figure;
            hold on;plot(BestSelfCSnorm(indLes([1:3 10:11])),BestLMANCSnorm(indLes([1:3 10:11])),'.','MarkerSize',25,'Color','r')
            plot(BestSelfCSnorm(indLes(4:9)),BestLMANCSnorm(indLes(4:9)),'.','MarkerSize',25,'Color','b')
            plot(BestSelfCSnorm(indLes(12:13)),BestLMANCSnorm(indLes(12:13)),'.','MarkerSize',25,'Color','g')
            hold on;plot(t,t,'-','Color','k')
            xlim([0 0.2]);ylim([0 0.2])
            plot([mean(BestSelfCSnorm(indLes))-std(BestSelfCSnorm(indLes))/sqrt(13);mean(BestSelfCSnorm(indLes))+std(BestSelfCSnorm(indLes))/sqrt(13)],[mean(BestLMANCSnorm(indLes));mean(BestLMANCSnorm(indLes))],'Color','k')
           plot([mean(BestSelfCSnorm(indLes));mean(BestSelfCSnorm(indLes))],[mean(BestLMANCSnorm(indLes))-std(BestLMANCSnorm(indLes))/sqrt(13);mean(BestLMANCSnorm(indLes))+std(BestLMANCSnorm(indLes))/sqrt(13)],'Color','k')
           
% Sigma calculation
    % see 820.m
    % or LOAD /cardinal/FiguresA/ConsensusData.mat
    
% Figure 9 - Zebra Finch Lesion Data
% /cardinal/FiguresA/FinalFigures/figure9data.mat
        clear crosscoAC
        clear crosscoINA
        clear ccAC
        clear ccINA
        clear mnccAC
        clear mnccINA
        clear mnccAC2
        clear mnccINA2
        for k=1:11
            i=k;
            notelength=PX(i,2)-PX(i,1);
            numms=floor(notelength/8);
            clear crosscoAC
            clear crosscoINA
            for ii=1:numms
                first=ii*8;
                middle=PX(i,1)+first;
                init=500-first+1;
                for j=1:notelength
                    ab=corrcoef(ZFresidsCTL(i).data(j+PX(i,1),:),ZFresidsCTL(i).data(middle,:));
                    crosscoAC(init+j,ii)=ab(2);
                    ai=corrcoef(ZFresidsINA(i).data(j+PX(i,1),:),ZFresidsINA(i).data(middle,:));
                    crosscoINA(init+j,ii)=ai(2);
                end
            end
            for i=1:size(crosscoAC,1)
                ind1=find(crosscoAC(i,:)>0);
                mnccAC(i,k)=mean(crosscoAC(i,ind1));
                ind2=find(crosscoINA(i,:)>0);
                mnccINA(i,k)=mean(crosscoINA(i,ind2));
            end
        end
        %% THIS TAKES THE MEAN
        for i=1:size(mnccAC,1)
            ind1=find(mnccAC(i,:)>0);
            mnccAC2(i)=mean(mnccAC(i,ind1).^2);
            stdAC2(i)=std(mnccAC(i,ind1).^2);
            ind2=find(mnccINA(i,:)>0);
            mnccINA2(i)=mean(mnccINA(i,ind1).^2);
            stdINA2(i)=std(mnccINA(i,ind1).^2);
        end
        sdAC=stdAC2/sqrt(11);
        sdINA=stdINA2/sqrt(11);
   %% THIS PLOTS IT
        t1=[1:1:length(mnccAC2)];
        figure;hold on;
        plot([t1;t1],[(mnccAC2-sdAC);(mnccAC2+sdAC)],'-','Color','r')
        plot([t1;t1],[(mnccINA2-sdINA);(mnccINA2+sdINA)],'-','Color','b')
        xlim([220 780]);ylim([0 1])
%%% FIGURES FOR LMAN only component
        for i=1:11
            CVratio(i)=1-mean(std(ZFresidsINA(i).data(PX(i,1):PX(i,2),:)'))/mean(std(ZFresidsCTL(i).data(PX(i,1):PX(i,2),:)'));
        end
        for i=1:size(mnccAC,2)
            weightedFF(:,i)=(mnccAC(:,i).^2-CVratio(i)*mnccINA(:,i).^2)/(1-CVratio(i));
            weightedFF(:,i)=weightedFF(:,i)/max(weightedFF(:,i));
        end
        for i=1:size(mnccAC,1)
            ind1=find(weightedFF(i,:)>0);
            weightedFF2(i)=mean(weightedFF(i,ind1));
        end
        weightedFF2=weightedFF2/max(weightedFF2);
        figure;plot(weightedFF2,'Linewidth',2)
        halfwidth=(500-399-min(find(weightedFF2(400:end)>0.5)))/4;
        % 13.215 (precisely)
        
        
        
%         for i=1:11
%             widthAC(i)=mean([mnccAC(430,i);mnccAC(270,i)]);
%             widthINA(i)=mean([mnccINA(430,i);mnccINA(270,i)]);
%         end
% figure;plot([1;2],[widthAC.^2;widthINA.^2],'k')
% hold on;plot(1,widthAC.^2,'.','MarkerSize',15,'Color','r')
% plot(2,widthINA.^2,'.','MarkerSize',15,'Color','b')
% xlim([0.5 2.5])
% ylim([0 1])
 

% FIGURE 9D - see correctforn.m
                  % Figure 9D
                %         t1=[1:1:size(mnBothAC,2)];
                %         figure;hold on;
                %         plot([t1;t1],[(mnBothINA-stdBothINA).^2;(mnBothINA+stdBothINA).^2],'-','Color','b')
                %         plot([t1;t1],[(mnBothAC-stdBothAC).^2;(mnBothAC+stdBothAC).^2],'-','Color','r')
                %         xlim([7 692]);ylim([0 1])
                
% FIGURE 10 - MSBanalyInfoTheory --- October 7, 2009
count=0;
for j=0.4:0.1:1
    count=count+1;
    for i=1:28

        if isequal(Predict(i).direction,'up')
            valu=60;
        else
            valu=40;
        end
        if i<23
            [pretime(count,i) posttime(count,i)]=ContingSimMSB(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).onset+400,:),valu,j);
        else
        	[pretime(count,i) posttime(count,i)]=ContingSimMSB(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),valu,j);
        end
    end
end

predicted=mean(posttime'-pretime')/8;
for i=1:22
    notewidth(i)=400;
end
for i=23:28
    notewidth(i)=Predict(i).offset-Predict(i).onset;
end
actual=mean(notewidth)/8;

[h,p]=ttest(posttime(5,:)-pretime(5,:),notewidth)
%   p=2.3e-8
for j=1:500
    predPrecise(j)=length(find(pretime(5,:)<-500+j))/28;
end
for j=501:1000
    predPrecise(j)=length(find(posttime(5,:)>j-500))/28;
end
figure;hold on;plot([-40:1/8:40],predPrecise(181:821))


% February 1, 2010

% This gives predictions at all points, -1 if no data
for i=1:22
    FFnorms(i,:)=ContingSimMSB2(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),80,Predict(i).direction);
end
%
for i=1:size(FFnorms,2)
    g=find(FFnorms(:,i)>-1);
    if isempty(g)
        FFnormmn(i)=0;
    else
        FFnormmn(i)=mean(FFnorms(g,i));
    end
end
%%%%

59823
for i=1:28
    FFnorms30(i,:)=ContingSimMSB2(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),80,Predict(i).direction,30,90);
end
%
for i=1:size(FFnorms30,2)
    g=find(FFnorms30(:,i)>-1);
    if isempty(g)
        FFnormmn30(i)=0;
    else
        FFnormmn30(i)=mean(FFnorms30(g,i));
    end
end
figure;plot(FFnormmn30,'b')
hold on;plot(FFnormmn20,'k')
hold on;plot(FFnormmn10,'r')
xlim([300 700])

%%%%%
% October 23, 2009
% Example with Targeting in the "middle"
% Experiment #30 - pk37bk19
% start at 290, go to 640

% This figure shows PRE/POST with overlayed contours
load 102309b.mat
figure;imagesc(t,f,log(avPre));syn
hold on;plot([0.080:1/8000:0.12],mean(pitchBaseline(769:1089,:)'),'.','Color','b')
xlim([0.08 0.12]);ylim([0 1e4])
figure; 
% pitchBaseline - all day 913
% pitchShifted - all day 921
diff=(mean(pitchpk37post(769:1089,:)')-mean(pitchpk37pre(769:1089,:)'))
figure;plot(diff/max(diff))
targ=(targeting+240)-32
hold on; plot([median(targ) median(targ)],[0 1])

%%%
% example 15
predMax=[zeros(1,135) ones(1,182) zeros(1,60)];
for i=1:length(predMax)-40
    predMaxSmooth(i)=mean(predMax(i:i+40));
end
%%% figure 4, diff pred diff
% figure;subplot(133);plot(Predict(24).LearnedNorm(203:483))
% hold on;plot(Predict(30).LearnedNorm(157:437))
% hold on;plot(CSs64(30,30).data(157:437)/0.016,'r')
% hold on;plot(CSs64(24,24).data(203:483)/0.01555,'r')
% hold on;plot([60 60],[0 1])
% xlim([20 220])
% ylim([0 1.1])
% 
% subplot(131);plot((Predict(24).ResidAC(223:423,1:31)),'b')

% 11.16.09
% SELF vs OTHER figure /mutual
figure;hold on;subplot(131);hold on
for i=1:31
    plot(Predict(24).ResidAC(243:423,i)-mean(Predict(24).ResidAC(243:423,1:31)')','b')
end
xlim([20 200]);ylim([-0.06 0.06])
subplot(132);
hold on;figure;hold on
for i=1:33
    plot(Predict(7).ResidAC(345:525,i)-mean(Predict(7).ResidAC(345:525,31:63)')','g')
end

xlim([20 200])
ylim([-0.06 0.06])

%%%%%%%%%5
%%%%%%%%%

for k=15
            i=k;
            notelength=200;
            numms=floor(notelength/4);
            clear crosscoAC
            clear crosscoINA
            for ii=1:numms
                first=ii*4;
                PX=860;
                middle=PX+first;
                init=300-first+1;
                for j=1:notelength
                    ab=corrcoef(Predict(i).ResidAC(j+PX,:),Predict(i).ResidAC(middle,:));
                    crosscoAC(init+j,ii)=ab(2);
                end
            end
            for i=1:size(crosscoAC,1)
                ind1=find(crosscoAC(i,:)>0);
                mn15(i)=mean(crosscoAC(i,ind1));
            end
end

figure;plot(mn24.^2);hold on;plot(mn15.^2,'r')

%%%%
BestSelfCSsyll64=[mean(BestSelfCSnorm64(1:2)) BestSelfCSnorm64(3) mean(BestSelfCSnorm64(4:6)) BestSelfCSnorm64(7) mean(BestSelfCSnorm64(8:9)) BestSelfCSnorm64(10:end)];
BestOtherCSsyll64=[mean(BestOtherCSnorm64(1:2)) BestOtherCSnorm64(3) mean(BestOtherCSnorm64(4:6)) BestOtherCSnorm64(7) mean(BestOtherCSnorm64(8:9)) BestOtherCSnorm64(10:end)];
figure;plot([1;2],[BestOtherCSsyll64;BestSelfCSsyll64],'k')
hold on;plot(1,BestOtherCSsyll64,'.','Markersize',15,'Color','b')
plot(2,BestSelfCSsyll64,'.','Markersize',15,'Color','r')
xlim([0.7 2.3])

% #7 and #24 - all this bullshit about #30 is ridiculous....
figure;hold on;
plot(CSs64(7,7).data(373:553)/0.0175)
plot(CSs64(24,24).data(243:423)/0.017/0.99/0.922)
plot(Predict(7).LearnedNorm(373:553),'r')
plot(Predict(24).LearnedNorm(243:423),'r')
ylim([0 1.05])
xlim([0 180])




%%%%%%%%%            figure;hold on;
            % Fig.4A - example
            figure;hold on
            subplot(121);hold on
            for i=1:20
            plot(Predict(7).ResidAC(Predict(7).onset:Predict(7).onset+350,i)-mean(Predict(7).ResidAC(Predict(7).onset:Predict(7).onset+350,1:20)')','r')
            end
            xlim([0 350]);ylim([-0.06 0.06])
            subplot(122);hold on
            for i=1:20
            plot(Predict(7).ResidINA(Predict(7).onset:Predict(7).onset+350,i)-mean(Predict(7).ResidINA(Predict(7).onset:Predict(7).onset+350,1:20)')','b')
            end
            xlim([0 350]);ylim([-0.06 0.06])
            
            % Add to figure 3E
            % hold on;plot([-25:50/400:25],FFnormmn10(300:700)/max(FFnormmn10),'Linewidth',3)
            
            % Add to figure 3C ---   355=500-[median(Predict(15).Targeting)-750]
            % hold on;plot(FFnorms10(15,355:750),'b')
            
% February 17, 2010            
59823
% First run figure 3G

% Then get the minimal one
                count=0;
                for j=0.4:0.1:1
                    count=count+1;
                    for i=1:28

                        if isequal(Predict(i).direction,'up')
                            valu=60;
                        else
                            valu=40;
                        end
                        if i<23
                            [pretime(count,i) posttime(count,i)]=ContingSimMSB(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).onset+400,:),valu,j);
                        else
                            [pretime(count,i) posttime(count,i)]=ContingSimMSB(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),valu,j);
                        end
                    end
                end

                %   p=2.3e-8
                for j=1:500
                    predPrecise(j)=length(find(pretime(5,:)<-500+j))/28;
                end
                for j=501:1000
                    predPrecise(j)=length(find(posttime(5,:)>j-500))/28;
                end
        hold on;plot([-40:1/8:40],predPrecise(181:821))
%%% Then factor in the peripheral constraints
                for kk=62
                for i=1:28
                    FFn(i,:)=ContingSimMSB2(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),80,Predict(i).direction,kk,90);
                end
                %
                for i=1:size(FFn,2)
                    g=find(FFn(:,i)>-1);
                    if isempty(g)
                        FFnmn(i)=0;
                    else
                        FFnmn(i)=mean(FFn(g,i));
                    end
                end
                distancefrom(kk)=sum(abs((mnabb(500:900)/max(mnabb(500:900)))-FFnmn(300:700)/max(FFnmn(300:700))));
                end
        hold on;plot(t(43:end)/8,FFnormmn60(1:1060)/max(FFnormmn60),'b')
        
        
% February 25, 2010
            figure;hold on;
            %             subplot(241);plot(Predict(7).ResidAC(Predict(7).onset+30:Predict(7).onset+350,70:90))
            %             xlim([0 320]);ylim([-0.06 0.06])
            mtar=mean(Predict(12).Targeting)
            mprc=prctile(Predict(12).ResidAC(round(mtar),:),70)

            indH=find(Predict(12).ResidAC(round(mtar),101:end)>mprc)
            indD=find(Predict(12).ResidAC(round(mtar),101:end)>-0.003 & Predict(12).ResidAC(round(mtar)-192,101:end)<0.003);
            subplot(211);hold on;plot(Predict(12).ResidAC(:,101:end),'k','Linewidth',2)
            hold on;plot(Predict(12).ResidAC(:,100+indH),'r','Linewidth',2)
            xlim([400 800])
            plot([400 800],[0 0],'k','Linewidth',4)
            subplot(212);hold on;plot(Predict(12).ResidAC(:,101:end),'k','Linewidth',2)
            hold on;plot(Predict(12).ResidAC(:,100+indD),'r','Linewidth',2)
            xlim([400 800])
            plot([400 800],[0 0],'k','Linewidth',4)
            

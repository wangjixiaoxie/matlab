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
      % AD7=Alldatasecondseven1sig(7);
        % createfigureA1A.m
        distance=112; % 30ms-16ms - distance between onset of spec and onset of contour
        initpt=250;
        finpt=900;
        % third harmonic
        meanBase=(mean(AD7.baselineAC'));
        meanShifted=(0.5*(mean(AD7.exp(17).selectedpitchcurves')+mean(AD7.exp(18).selectedpitchcurves')));
        createfigureA1A(log(OR92heatmap.baseline))
        xlim([-0.01 0.11]);ylim([0 10000])
        % ADJUST limits on colormap to 6/15
        hold on;plot(t1(initpt-distance:finpt-distance),meanBase(initpt:finpt),'LineWidth',2)
        %%
        createfigureA1A(log(OR92heatmap.shifted))
        xlim([-0.01 0.11]);ylim([0 10000])
        % ADJUST limits on colormap to 6/15
        hold on;plot(t1(initpt-distance:finpt-distance),meanShifted(initpt:finpt),'LineWidth',2)
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
            ylim([0 1.05])
            xlim([-40 40])


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
        plot([middle],tops,'.','MarkerSize','k')
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

% Figure 4 - inter-syllable comparison
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
% Figure 5 -                 
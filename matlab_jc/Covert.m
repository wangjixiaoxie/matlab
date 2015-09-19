% Figure 1B - Block of LMAN-->RA (expression of adaptation)


% Figure 3A - changes are specific - sub-syllabic localization
figure;hold on;
count=0;
for i=[1 7 10:11]
    window(i,:)=round([median(Experiment(i).TargetingWN)-64 median(Experiment(i).TargetingWN)]);
    window(7,:)=[147 170];

    count=count+1;
    subplot(2,2,count);hold on;
    clear timedist
    for j=2:length(Experiment(i).timeACpost)
        timedist(j-1)=Experiment(i).timeACpost(j)-Experiment(i).timeACpost(j-1);
    end
    clear ind
    ind=find(timedist>8)+1;
    if i==7
        Adchange(i).data=(mean(Experiment(i).pitchACpost(:,ind(1):end)')-mean(Experiment(i).pitchACpre(:,50:end)'));
        plot(Adchange(i).data)
        VarRecov(i)=mean(std(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):end)'))/mean(std(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,end-20:end)'))
        xlim([135 420])
    else if i==10
            Adchange(i).data=(mean(Experiment(i).pitchACpost(:,ind(2):end)')-mean(Experiment(i).pitchACpre(:,50:end)'));
            plot(Adchange(i).data)
            VarRecov(i)=mean(std(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,ind(2):end)'))/mean(std(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,end-50:end)'))
            xlim([150 450])
        else if i==11
               Adchange(i).data=(mean(Experiment(i).pitchACpost(:,ind(2):ind(3)-1)')-mean(Experiment(i).pitchACpre(:,50:end)'));
               plot(Adchange(i).data)
               VarRecov(i)=mean(std(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,ind(2):ind(3)-1)'))/mean(std(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,end-50:end)'))
               xlim([150 450])
            else
               Adchange(i).data=(mean(Experiment(i).pitchACpost(:,ind(1):ind(2)-1)')-mean(Experiment(i).pitchACpre(:,50:end)'));
               VarRecov(i)=mean(std(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,ind(1):ind(2)-1)'))/mean(std(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,end-50:end)'))
               plot(Adchange(i).data)
               xlim([130 330]);ylim([0 110])
            end
        end
    end
    [b,a]=hist(Experiment(i).TargetingWN(find(Experiment(i).TargetingWN<600))-32,10);
    stairs(a,(b/sum(b))*150,'r')
end
% Figure 1A - Reversible block of LMAN-->RA (variability)
        for i=[1:3 5:11]
            window(i,:)=round([median(Experiment(i).TargetingWN)-64 median(Experiment(i).TargetingWN)]);
            VarReduction(i)=mean(std(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,end-50:end)'))/mean(std(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,end-50:end)'));
            % first 30 notes is just onset of APV effect
            if i==7
                VarReduction(i)=mean(std(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,end-20:end)'))/mean(std(Experiment(10).pitchACpre(Experiment(i).on:Experiment(i).off,end-20:end)'));
            end
            % not the first experiment
            if i==8 || i==9
                VarReduction(i)=mean(std(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,end-50:end)'))/mean(std(Experiment(7).pitchACpre(Experiment(i).on:Experiment(i).off,end-50:end)'));
            end
            if i==11
                VarReduction(i)=mean(std(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,end-50:end)'))/mean(std(Experiment(10).pitchACpre(Experiment(i).on:Experiment(i).off,end-50:end)'));
            end
        end
        figure;plot([1;2;3],[ones(1,4);VarReduction([1 7 10:11]);VarRecov([1 7 10 11])],'Color','k')
        hold on;plot(1,1,'+','MarkerSize',20)
        plot(2,VarReduction([1 7 10:11]),'+','MarkerSize',20,'Color','b')
        plot(3,VarRecov([1 7 10:11]),'+','MarkerSize',20,'Color','b')
        
        
% Figure 1 - example - experiment 11
i=7;
window(i,:)=round([median(Experiment(i).TargetingWN)-64 median(Experiment(i).TargetingWN)]);
window(7,:)=[147 170];
figure;hold on;
plot(Experiment(i).timeACpre,mean(Experiment(i).pitchACpre(window(i,1):window(i,2),:)),'.','Markersize',6,'Color','b')
plot(Experiment(i).timeAPV,mean(Experiment(i).pitchAPV(window(i,1):window(i,2),:)),'.','Markersize',6,'Color','g')
plot(Experiment(i).timeAPVwn,mean(Experiment(i).pitchAPVwn(window(i,1):window(i,2),:)),'.','Markersize',6,'Color','r')
plot(Experiment(i).timeACpost,mean(Experiment(i).pitchACpost(window(i,1):window(i,2),:)),'.','Markersize',6,'Color','b')

mean(std(Experiment(i).pitchAPV(window(i,1):window(i,2),end-50:end)'))/mean(std(Experiment(i).pitchACpre(window(i,1):window(i,2),:)'))
% for i=11
hold on;plot([5645 5825],[ mean(mean(Experiment(i).pitchACpre(window(i,1):window(i,2),:))) mean(mean(Experiment(i).pitchACpre(window(i,1):window(i,2),:)))],'k')
% reduction = 47% --- WOW!
% bad example b/c no FF reduction in this note (although there is a ton in
% the next note!!


% Figure 2
 load /cardinal6/covertAnaly/Summary91009.mat

    % PRE1 - mean of data for ~2 days prior to shift
    % PRE2 - mean of pitch across entire targ region (unless it is noisy (i.e. WN), in which case pitch only in the non-noisy part of the targ region)
    % APV1 - mean pitch for first half of notes with APVwn (large enough sample
    %           size, and if there is a real change, then any adaptive change
    %           would cause the first half to be less adaptive than the second
    %           half --- by eliminating any exps where this happens (red
    %           traces), we can ensure we only look at exps with no adaptive
    %           change.  Caveat is that there occassionally are adaptive
    %           changes near the beginning or end, but in ALL the black traces,
    %           they are counteracted by maladaptive changes elsewhere.
    % APV2 - mean pitch for second half of notes with APVwn
    % POST1 - mean pitch across entire DAY n+1 (where day n=WNoff, ACSF on)
    % POST 2 - mean pitch at the end of the experiment (missing n=2 because did
    %           incremental afterwards OR the bird died)

            clear SPlot
            clear Plotter
            clear PlotterNorm
            for i=[1:7 10 12:13]
                    window(i,:)=round([median(Experiment(i).TargetingWN)-64 median(Experiment(i).TargetingWN)]);
                    window(7,:)=[147 170]; % because it gets noisy afterwards - WN
                    window(4,:)=[200 270];% because of WN
                    SPlot(i).PRE1=mean(mean(Experiment(i).pitchACpre(window(i,1):window(i,2),:)'));
                    SPlot(i).PRE2=mean(mean(Experiment(i).pitchACpre(window(i,1):window(i,2),end-50:end)'));
                    llength=round(size(Experiment(i).pitchAPVwn,2)/2);
                    SPlot(i).APV1=mean(mean(Experiment(i).pitchAPV(window(i,1):window(i,2),:)'));
                    SPlot(i).APV2=mean(mean(Experiment(i).pitchAPVwn(window(i,1):window(i,2),end-llength:end)'));
                    %SPlot(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),end-50:end)'));
                    clear timedist
                    for j=2:length(Experiment(i).timeACpost)
                        timedist(j-1)=Experiment(i).timeACpost(j)-Experiment(i).timeACpost(j-1);
                    end
                    %figure;plot(timedist)
                    %%% entire next day
                    clear ind
                    ind=find(timedist>8)+1;

                        SPlot(i).POST1=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),1:ind(1))'));
                        if size(ind,2)==1
                            SPlot(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):end)'));
                        else
                            if size(ind,2)==2
                                SPlot(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):ind(2))'));
                                SPlot(i).POST3=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(2):end)'));
                            else
                                SPlot(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):ind(2))'));
                                SPlot(i).POST3=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(2):ind(3))'));
                            end
                        end

            end

            znorm=[21 21 0 0 0 0 30 0 0 26.3 0 12 12];
            goodinas=[1:7 10 12:13];
                for i=goodinas
                    if isequal(Experiment(i).DIR,'up')
                        coef=1;
                    else
                        coef=-1;
                    end
                    Plotter(i,1)=coef*((SPlot(i).PRE1-SPlot(i).PRE2));
                    Plotter(i,2)=0;
                    Plotter(i,3)=coef*((SPlot(i).APV1-SPlot(i).PRE2));
                    Plotter(i,4)=coef*((SPlot(i).APV2-SPlot(i).PRE2));
                    Plotter(i,5)=coef*((SPlot(i).POST1-SPlot(i).PRE2)); 
                    if ~isempty(SPlot(i).POST2)
                        Plotter(i,6)=coef*((SPlot(i).POST2-SPlot(i).PRE2));
                    end
                    if ~isempty(SPlot(i).POST3)
                        Plotter(i,7)=coef*((SPlot(i).POST3-SPlot(i).PRE2));
                    end
                    PlotterNorm(i,1)=coef*((SPlot(i).PRE1-SPlot(i).PRE2))*(1/znorm(i));
                    PlotterNorm(i,2)=0;
                    PlotterNorm(i,3)=coef*((SPlot(i).APV1-SPlot(i).PRE2))*(1/znorm(i));
                    PlotterNorm(i,4)=coef*((SPlot(i).APV2-SPlot(i).PRE2))*(1/znorm(i));
                    PlotterNorm(i,5)=coef*((SPlot(i).POST1-SPlot(i).PRE2))*(1/znorm(i)); 
                    if ~isempty(SPlot(i).POST2)
                        PlotterNorm(i,6)=coef*((SPlot(i).POST2-SPlot(i).PRE2))*(1/znorm(i));
                    end
                    if ~isempty(SPlot(i).POST3)
                        PlotterNorm(i,7)=coef*((SPlot(i).POST3-SPlot(i).PRE2))*(1/znorm(i));
                    end
                end


                  
                
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%



        clear SPlotR
        clear PlotterR
        clear PlotterNormR
            for i=[1:7 10 12:13]
                    windowCTL=[Experiment(i).onCTL Experiment(i).offCTL];
                    window(i,:)=round([median(Experiment(i).TargetingWN)-64 median(Experiment(i).TargetingWN)]);
                    window(7,:)=[147 170]; % because it gets noisy afterwards - WN
                    window(4,:)=[200 270];% because of WN
                    SPlotR(i).PRE1=mean(mean(Experiment(i).pitchACpre(window(i,1):window(i,2),:)'))-mean(mean(Experiment(i).pitchACpreCTL(windowCTL,:)'));
                    SPlotR(i).PRE2=mean(mean(Experiment(i).pitchACpre(window(i,1):window(i,2),end-100:end)'))-mean(mean(Experiment(i).pitchACpreCTL(windowCTL,end-100:end)'));
                    llength=round(size(Experiment(i).pitchAPVwn,2)/2);
                    SPlotR(i).APV1=mean(mean(Experiment(i).pitchAPV(window(i,1):window(i,2),:)'))-mean(mean(Experiment(i).pitchAPVCTL(windowCTL,:)'));
                    SPlotR(i).APV2=mean(mean(Experiment(i).pitchAPVwn(window(i,1):window(i,2),end-llength:end)'))-mean(mean(Experiment(i).pitchAPVwnCTL(windowCTL,end-llength:end)'));
                    %SPlot(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),end-50:end)'));
                    clear timedist
                    for j=2:length(Experiment(i).timeACpost)
                        timedist(j-1)=Experiment(i).timeACpost(j)-Experiment(i).timeACpost(j-1);
                    end
                    %figure;plot(timedist)
                    %%% entire next day
                    clear ind
                    ind=find(timedist>8)+1;

                        SPlotR(i).POST1=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),1:ind(1))'))-mean(mean(Experiment(i).pitchACpostCTL(windowCTL,1:ind(1))'));
                        if size(ind,2)==1
                            SPlotR(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):end)'))-mean(mean(Experiment(i).pitchACpostCTL(windowCTL,ind(1):end)'));
                        else
                            if size(ind,2)==2
                                SPlotR(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):ind(2))'))-mean(mean(Experiment(i).pitchACpostCTL(windowCTL,ind(1):ind(2))'));
                                SPlotR(i).POST3=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(2):end)'))-mean(mean(Experiment(i).pitchACpostCTL(windowCTL,ind(2):end)'));
                            else
                                SPlotR(i).POST2=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):ind(2))'))-mean(mean(Experiment(i).pitchACpostCTL(windowCTL,ind(1):ind(2))'));
                                SPlotR(i).POST3=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(2):ind(3))'))-mean(mean(Experiment(i).pitchACpostCTL(windowCTL,ind(2):ind(3))'));
                            end
                        end

            end
            
            znorm=[18 18 0 0 0 0 15 0 0 19.8 0 6 6];
                goodinas=[1:7 10 12:13];
                for i=goodinas
                    if isequal(Experiment(i).DIR,'up')
                        coef=1;
                    else
                        coef=-1;
                    end
                    PlotterNormR(i,1)=coef*((SPlotR(i).PRE1-SPlotR(i).PRE2))*(1/znorm(i));
                    PlotterNormR(i,2)=0;
                    PlotterNormR(i,3)=coef*((SPlotR(i).APV1-SPlotR(i).PRE2))*(1/znorm(i));
                    PlotterNormR(i,4)=coef*((SPlotR(i).APV2-SPlotR(i).PRE2))*(1/znorm(i));
                    PlotterNormR(i,5)=coef*((SPlotR(i).POST1-SPlotR(i).PRE2))*(1/znorm(i)); 
                    if ~isempty(SPlot(i).POST2)
                        PlotterNormR(i,6)=coef*((SPlotR(i).POST2-SPlotR(i).PRE2))*(1/znorm(i));
                    end
                    if ~isempty(SPlotR(i).POST3)
                        PlotterNormR(i,7)=coef*((SPlotR(i).POST3-SPlotR(i).PRE2))*(1/znorm(i));
                    end
                    PlotterR(i,1)=coef*((SPlotR(i).PRE1-SPlotR(i).PRE2));
                    PlotterR(i,2)=0;
                    PlotterR(i,3)=coef*((SPlotR(i).APV1-SPlotR(i).PRE2));
                    PlotterR(i,4)=coef*((SPlotR(i).APV2-SPlotR(i).PRE2));
                    PlotterR(i,5)=coef*((SPlotR(i).POST1-SPlotR(i).PRE2)); 
                    if ~isempty(SPlotR(i).POST2)
                        PlotterR(i,6)=coef*((SPlotR(i).POST2-SPlotR(i).PRE2));
                    end
                     if ~isempty(SPlotR(i).POST3)
                         PlotterR(i,7)=coef*((SPlotR(i).POST3-SPlotR(i).PRE2));
                     end
                end

                
                figure;hold on;subplot(223);hold on;
                plot([0 6],[0 0],'r','Linewidth',2)
                figure;plot([1:4 5],[PlotterR(goodinas,[1:4 6])'],'Linewidth',2,'Color','k')
                plot(1,PlotterR([goodinas],1),'.','MarkerSize',20,'Color','b')
                plot(2,PlotterR([goodinas],2),'.','MarkerSize',20,'Color','b')
                plot(3,PlotterR([goodinas],3),'.','MarkerSize',20,'Color','b')
                plot(4,PlotterR([goodinas],4),'.','MarkerSize',20,'Color','b')
                plot(5,PlotterR([goodinas],6),'.','MarkerSize',20,'Color','b')
                ind7=find(PlotterR(:,7)~=0);
                plot([1:6],[PlotterR(ind7,[1:4 6:7])'],'Linewidth',2,'Color','k')
                xlim([0.8 6.2])
                subplot(224);hold on;
                plot([0 6],[0 0],'r','Linewidth',2)
                plot([1:4 5],[PlotterNormR(goodinas,[1:4 6])'],'Linewidth',2,'Color','k')
                plot(1,PlotterNormR([goodinas],1),'.','MarkerSize',20,'Color','b')
                plot(2,PlotterNormR([goodinas],2),'.','MarkerSize',20,'Color','b')
                plot(3,PlotterNormR([goodinas],3),'.','MarkerSize',20,'Color','b')
                plot(4,PlotterNormR([goodinas],4),'.','MarkerSize',20,'Color','b')
                plot(5,PlotterNormR([goodinas],6),'.','MarkerSize',20,'Color','b')
                ind7=find(PlotterNormR(:,7)~=0);
                plot([1:6],[PlotterNormR(ind7,[1:4 6:7])'],'Linewidth',2,'Color','k')

                xlim([0.8 6.2])


              figure; subplot(221);hold on;
                plot([0 6],[0 0],'r','Linewidth',2)
                plot([1:5],[Plotter(goodinas,[1:4 6])'],'Linewidth',2,'Color','k')
                plot(1,Plotter([goodinas],1),'.','MarkerSize',20,'Color','b')
                plot(2,Plotter([goodinas],2),'.','MarkerSize',20,'Color','b')
                plot(3,Plotter([goodinas],3),'.','MarkerSize',20,'Color','b')
                plot(4,Plotter([goodinas],4),'.','MarkerSize',20,'Color','b')
                plot(5,Plotter([goodinas],6),'.','MarkerSize',20,'Color','b')
                ind7=find(Plotter(goodinas,7)~=0);
                plot([1:6],[Plotter(goodinas(ind7),[1:4 6:7])'],'Linewidth',2,'Color','k')
                xlim([0.8 6.2])
                subplot(222);hold on;
                plot([0 6],[0 0],'r','Linewidth',2)
                plot([1:5],[PlotterNorm(goodinas,[1:4 6])'],'Linewidth',2,'Color','k')
                plot(1,PlotterNorm([goodinas],1),'.','MarkerSize',20,'Color','b')
                plot(2,PlotterNorm([goodinas],2),'.','MarkerSize',20,'Color','b')
                plot(3,PlotterNorm([goodinas],3),'.','MarkerSize',20,'Color','b')
                plot(4,PlotterNorm([goodinas],4),'.','MarkerSize',20,'Color','b')
                plot(5,PlotterNorm([goodinas],6),'.','MarkerSize',20,'Color','b')
                ind7=find(PlotterNorm(goodinas,7)~=0);
                plot([1:6],[PlotterNorm(goodinas(ind7),[1:4 6:7])'],'Linewidth',2,'Color','k')
                xlim([0.8 6.2])
                





%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                

                
% Fairly robust
     % AC1 - pitchvals in pvs 12hrs (robust to other values e.g. 4)
     % APV1 - pitchvals w/ APV on but before wn on (last 20 notes - near end)
     % APV2 - pitchvals over the last 4hrs of songs with APV and wn on for
            % 24hr shifts (~14hrs of learning) and last quarter of notes for
            % faster shfits (thus ~ last 25% of learning for all exps)
     % AC2 - pitchvals in first 4hrs
            for i=1:6
                indA1(i)=min(find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-12));

                ons=Experiment(i).on;
                offs=Experiment(i).off;
                AC1(i)=median(median(Experiment(i).pitchACpre(ons:offs,indA1(i):end)'));
                allAPV=Experiment(i).pitchAPVwn;
                llength=size(allAPV,2);
                indB1(i)=min(find(Experiment(i).timeAPVwn>Experiment(i).timeAPVwn(end)-4));
                APV1(i)=median(median(Experiment(i).pitchAPV(ons:offs,end-20:end)'));
                APV2(i)=median(median(allAPV(ons:offs,indB1(i):end)'));
              % first 4 hrs after
                ind1(i)=min(find(Experiment(i).timeACpost>Experiment(i).timeAPVwn(end)+0));
                ind2(i)=min(find(Experiment(i).timeACpost>Experiment(i).timeAPVwn(end)+4));
                ind1(4)=37;
                ind2(4)=56;
                AC2(i)=median(median(Experiment(i).pitchACpost(ons:offs,ind1(i):ind2(i))'));
            end
            for i=7:8
                indA1(i)=min(find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-12));

                ons=Experiment(i).on;
                offs=Experiment(i).off;
                AC1(i)=median(median(Experiment(i).pitchACpre(ons:offs,indA1(i):end)'));
                allAPV=Experiment(i).pitchAPVwn;
                llength=size(allAPV,2);
                indB1(i)=min(find(Experiment(i).timeAPVwn>Experiment(i).timeAPVwn(end)-4));
                APV1(i)=median(median(Experiment(i).pitchAPV(ons:offs,end-20:end)'));
                APV2(i)=median(median(allAPV(ons:offs,end-llength/4:end)'));
              % first 4 hrs after
                ind1(i)=min(find(Experiment(i).timeACpost>Experiment(i).timeAPVwn(end)+0));
                ind2(i)=min(find(Experiment(i).timeACpost>Experiment(i).timeAPVwn(end)+4));
                ind1(4)=37;
                ind2(4)=56;
                AC2(i)=median(median(Experiment(i).pitchACpost(ons:offs,ind1(i):ind2(i))'));
            end
            for i=9
                indA1(i)=min(find(Experiment(i).timeACpre>Experiment(i).timeACpre(end)-12));
                ons=Experiment(i).on;
                offs=Experiment(i).off;
                AC1(i)=median(median(Experiment(i).pitchACpre(ons:offs,indA1(i):end)'));
                allAPV=Experiment(i).pitchAPVwn;
                llength=size(allAPV,2);
                APV1(i)=median(median(Experiment(i).pitchAPV(ons:offs,end-20:end)'));
                APV2(i)=median(median(allAPV(ons:offs,end-llength/4:end)'));
              % first 1-3 hrs after
                ind1(i)=min(find(Experiment(i).timeACpost>Experiment(i).timeAPVwn(end)+0));
                %ind2(i)=min(find(Experiment(i).timeACpost>Experiment(i).timeAPVwn(end)+4));
                ind1(4)=37;
                ind2(4)=56;
                AC2(i)=median(median(Experiment(i).pitchACpost(ons:offs,ind1(i):end)'));
            end

            AC1(3)=median(median(Experiment(3).pitchACpre(190:250,:)'));
            AC2(3)=median(median(Experiment(3).pitchACpost(190:250,ind1(3):ind2(3))'));
            for i=1:9
                if isequal(Experiment(i).DIR,'down')
                    fac=-1;
                else
                    fac=1;
                end
                % How far does the AFP push pitch in the adaptive direction before WN
                AFPpre(i)=fac*(AC1(i)-APV1(i));
                % How far does the AFP push pitch in the adaptive direction after WN
                AFPpost(i)=fac*(AC2(i)-APV2(i));
            end
            [h,p]=ttest(AFPpost([1:2 5:9])-AFPpre([1:2 5:9]))
            figure;hold on;
            subplot(122);hold on;plot([1 6 7 8 9],AFPpost([1 6 7 8 9])-AFPpre([1 6 7 8 9]),'.','MarkerSize',15)
            plot([2 5],AFPpost([2 5])-AFPpre([2 5]),'.','MarkerSize',15,'Color','r')
            hold on;plot([0 10],[0 0],'k')
            %ylim([-20 100])
            subplot(121);hold on;plot([AFPpre([1 2 5 6 7 8 9]);AFPpost([1 2 5 6 7 8 9])],'Color','k')
            plot(1,AFPpre([1 6 7 8 9]),'.','MarkerSize',15,'Color','b')
            plot(2,AFPpost([1 6 7 8 9]),'.','MarkerSize',15,'Color','b')
            plot(1,AFPpre([2 5]),'.','MarkerSize',15,'Color','r')
            plot(2,AFPpost([2 5]),'.','MarkerSize',15,'Color','r')
            xlim([0 3])

for i=1:9
    ons=Experiment(i).on-50;
    offs=Experiment(i).off+50;

    CVapv(i)=min(std(Experiment(i).pitchAPVwn(ons:offs,1:50)'));
    CVac(i)=min(std(Experiment(i).pitchACpre(ons:offs,:)'));
end
    
i=2
figure;plot(Experiment(i).timeACpre,median(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,:)),'*','Color','b')
hold on;plot(Experiment(i).timeAPV,median(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,:)),'*','Color','r')
hold on;plot(Experiment(i).timeAPVwn,median(Experiment(i).pitchAPVwn(Experiment(i).on:Experiment(i).off,:)),'*','Color','g')
hold on;plot(Experiment(i).timeACpost,median(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,:)),'*','Color','b')
hold on;plot(4978,median(median(Experiment(i).pitchACpre(Experiment(i).on:Experiment(i).off,:))),'+','MarkerSize',25)
hold on;plot(4986,median(median(Experiment(i).pitchACpost(Experiment(i).on:Experiment(i).off,:))),'+','MarkerSize',25,'Color','k')
hold on;plot(4979.5,median(median(Experiment(i).pitchAPV(Experiment(i).on:Experiment(i).off,:))),'+','MarkerSize',25,'Color','k')
hold on;plot(4984,median(median(Experiment(i).pitchAPVwn(Experiment(i).on:Experiment(i).off,70:124))),'+','MarkerSize',25,'Color','k')




                    if i==4
                        ind=find(2830<Experiment(4).timeACpost<2850);
                    else if i==5
                            ind=find(2910<Experiment(5).timeACpost<2930);
                        end
                    end
                    % This chooses the POST day when variability recovers - make this less ad hoc?   
                    if i==7
                        SPlot(7).POST1=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):end)'));
                    else if i==10
                            SPlot(10).POST1=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(2):end)'));
                    else if i==11
                            SPlot(11).POST1=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(2):ind(3)-1)'));
                        else
                            SPlot(i).POST1=mean(mean(Experiment(i).pitchACpost(window(i,1):window(i,2),ind(1):ind(2)-1)'));
                        end
                        end
                    end

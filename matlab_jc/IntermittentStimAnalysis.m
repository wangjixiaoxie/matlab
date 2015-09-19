%%%%%%%%%%%%%%
% "Dynamics of consolidation"
%%%%%%%%%%%%%%%
% Pretty straightforward, load the relevant file, daynight.m does the
% calculations very quickly, then plot the running avg'd data and linear
% fits for each day.
%%%%%%%%%%%%%%
% Things to keep in mind:
% 1. quality of reversion, since poor quality inactivations would be confounded by
% dynamics in the AFP component of learning (as opposed to dynamics of
% consolidation).
% 2. baseline dynamics
%%%%%%%%%%%%%%

% April 12,2011
load /bulbul1/LMAN_microstimulation/IntermittentAnaly041211.mat
IntExp(1).dir=1;IntExp(2).dir=1;IntExp(3).dir=1;IntExp(4).dir=-1;
for i=1:4
    % 1=all, 2=day, 3=night
    Learn{i}(1,:)=IntExp(i).dir*IntExp(i).deltadayFBlearning(1:length(IntExp(i).deltanightFBlearning))+IntExp(i).deltanightFBlearning;
    Learn{i}(2,:)=IntExp(i).dir*IntExp(i).deltadayFBlearning(1:length(IntExp(i).deltanightFBlearning));
    Learn{i}(3,:)=IntExp(i).dir*IntExp(i).deltanightFBlearning;    
end
% Experiment 1 is kind of crappy
figure;hold on;
for i=1:4
    subplot(2,2,i);hold on;
    plot(Learn{i}(1,:),Learn{i}(2,:),'r.','Markersize',15) % daytime
    plot(Learn{i}(1,:),Learn{i}(3,:),'k.','Markersize',15) % nighttime
    plot([-80 60],[-80 60],'k');plot([-80 60],[0 0],'k');plot([0 0],[-80 60],'k')
end

all=([Learn{2}(1,:) Learn{3}(1,:) Learn{4}(1,:)])
day=([Learn{2}(2,:) Learn{3}(2,:) Learn{4}(2,:)])
night=([Learn{2}(3,:) Learn{3}(3,:) Learn{4}(3,:)])
corrcoef(all,day)
corrcoef(day,night)
% Analysis of daynight
% r30g38 (December 2010) - decent experiment
                load /cardinal2/StimBirds/r30g38/experimentB.mat
              ExpIndex=1;
                IntExp(ExpIndex).bird='r30g38';    
                IntExp(ExpIndex).window=[300:340];
                IntExp(ExpIndex).FFdataFBlearning=[mean(experiment(10).contours(IntExp(ExpIndex).window,experiment(10).crfbind)) mean(experiment(11).contours(IntExp(ExpIndex).window,experiment(11).crfbind))];
                IntExp(ExpIndex).fvalsFBlearning=[experiment(10).fv(experiment(10).crfbind) experiment(11).fv(experiment(11).crfbind)];
                [IntExp(ExpIndex).deltadayFBlearning,IntExp(ExpIndex).deltanightFBlearning,IntExp(ExpIndex).AllDaysFBlearning,IntExp(ExpIndex).FFmorningfitFBlearning,IntExp(ExpIndex).FFeveningfitFBlearning]=daynight...
                    (IntExp(ExpIndex).FFdataFBlearning,IntExp(ExpIndex).fvalsFBlearning);


                IntExp(ExpIndex).FFdataCTlearning=[mean(experiment(10).contours(IntExp(ExpIndex).window,experiment(10).crctind)) mean(experiment(11).contours(IntExp(ExpIndex).window,experiment(11).crctind))];
                IntExp(ExpIndex).fvalsCTlearning=[experiment(10).fv(experiment(10).crctind) experiment(11).fv(experiment(11).crctind)];
                [IntExp(ExpIndex).deltadayCTlearning,IntExp(ExpIndex).deltanightCTlearning,IntExp(ExpIndex).AllDaysCTlearning,IntExp(ExpIndex).FFmorningfitCTlearning,IntExp(ExpIndex).FFeveningfitCTlearning]=daynight...
                    (IntExp(ExpIndex).FFdataCTlearning,IntExp(ExpIndex).fvalsCTlearning);    



                IntExp(ExpIndex).FFdataFBbaseline=mean(experiment(9).contours(IntExp(ExpIndex).window,experiment(9).crfbind(272:end)));
                IntExp(ExpIndex).fvalsFBbaseline=experiment(9).fv(experiment(9).crfbind(272:end));
                [IntExp(ExpIndex).deltadayFBbaseline,IntExp(ExpIndex).deltanightFBbaseline,IntExp(ExpIndex).AllDaysFBbaseline,IntExp(ExpIndex).FFmorningfitFBbaseline,IntExp(ExpIndex).FFeveningfitFBbaseline]=daynight...
                    (IntExp(ExpIndex).FFdataFBbaseline,IntExp(ExpIndex).fvalsFBbaseline);


                IntExp(ExpIndex).FFdataCTbaseline=mean(experiment(9).contours(IntExp(ExpIndex).window,experiment(9).crctind(1071:end)));
                IntExp(ExpIndex).fvalsCTbaseline=experiment(9).fv(experiment(9).crctind(1071:end));
                [IntExp(ExpIndex).deltadayCTbaseline,IntExp(ExpIndex).deltanightCTbaseline,IntExp(ExpIndex).AllDaysCTbaseline,IntExp(ExpIndex).FFmorningfitCTbaseline,IntExp(ExpIndex).FFeveningfitCTbaseline]=daynight...
                    (IntExp(ExpIndex).FFdataCTbaseline,IntExp(ExpIndex).fvalsCTbaseline);    

            % What is the reversion relative to baseline effect???
                        FBB=[];
                        for i=1:length(IntExp(ExpIndex).AllDaysFBbaseline)
                            FBB=[FBB IntExp(ExpIndex).AllDaysFBbaseline(i).FFdata];
                        end
                        CTB=[];
                        for i=1:length(IntExp(ExpIndex).AllDaysCTbaseline)
                            CTB=[CTB IntExp(ExpIndex).AllDaysCTbaseline(i).FFdata];
                        end
                        FBbaseFF=median(FBB);
                        CTbaseFF=median(CTB);
                        FBLearning=[];
                        CTLearning=[];
                        for i=1:length(IntExp(ExpIndex).AllDaysFBlearning)
                            FBLearning(i)=median(IntExp(ExpIndex).AllDaysFBlearning(i).FFdata);
                            CTLearning(i)=median(IntExp(ExpIndex).AllDaysCTlearning(i).FFdata);
                        end
                        FBLearning-median(FBB)
                        CTLearning-median(CTB)
                        1-(FBLearning-median(FBB))./(CTLearning-median(CTB)) % Reversion %age (0.0 is no reversion, 1.0 is complete reversion)


                %%%%%%%%%%%%%%%%%%%%%%%%%%5                      
                        figure;hold on;
                    for i=[9:11]% 14]%2:5%5:9%17:length(experiment)%[9:10 14:15]%length(experiment)
                        wide=20;
                        window=[300:340];
                        if size(experiment(i).crctind>wide)
                            [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                            timedata=b;
                            FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                            lastpoints=find(diff(timedata)>10);
                            startpoints=[1 lastpoints+1];
                            endpoints=[lastpoints length(timedata)];
                            for j=1:length(endpoints)
                                if endpoints(j)>startpoints(j)+wide-1
                                    if i==5 | i==6
                                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','b')
                                    else
                                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','b')
                                    end
                                end

                            end
                        end
                        if size(experiment(i).crfbind>wide)
                            [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                            timedata=b;
                            FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                            lastpoints=find(diff(timedata)>10);
                            startpoints=[1 lastpoints+1];
                            endpoints=[lastpoints length(timedata)];
                            for j=1:length(endpoints)
                                if endpoints(j)>startpoints(j)+wide-1
                                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','r')
                                end
                            end
                        end
                    end
                 %%%%%%%%%%%%%%%%
                 for i=1:length(IntExp(ExpIndex).FFmorningfitFBlearning)
                     t1=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(1));
                     t2=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(end));
                     plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBlearning(i);IntExp(ExpIndex).FFeveningfitFBlearning(i)],'k-','Linewidth',4)
                 end
                 


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% bk48w74 - very good experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    load /cardinal3/StimBirds/bk48w74/experiments.mat
                  ExpIndex=2;
                    IntExp(ExpIndex).bird='bk48w74';    
                    IntExp(ExpIndex).window=[300:350];
                    IntExp(ExpIndex).FFdataFBlearning=[mean(experiment(2).contours(IntExp(ExpIndex).window,experiment(2).crfbind)) mean(experiment(3).contours(IntExp(ExpIndex).window,experiment(3).crfbind))];
                    IntExp(ExpIndex).fvalsFBlearning=[experiment(2).fv(experiment(2).crfbind) experiment(3).fv(experiment(3).crfbind)];
                    [IntExp(ExpIndex).deltadayFBlearning,IntExp(ExpIndex).deltanightFBlearning,IntExp(ExpIndex).AllDaysFBlearning,IntExp(ExpIndex).FFmorningfitFBlearning,IntExp(ExpIndex).FFeveningfitFBlearning]=daynight...
                        (IntExp(ExpIndex).FFdataFBlearning,IntExp(ExpIndex).fvalsFBlearning);

                    IntExp(ExpIndex).FFdataCTlearning=[mean(experiment(2).contours(IntExp(ExpIndex).window,experiment(2).crctind)) mean(experiment(3).contours(IntExp(ExpIndex).window,experiment(3).crctind))];
                    IntExp(ExpIndex).fvalsCTlearning=[experiment(2).fv(experiment(2).crctind) experiment(3).fv(experiment(3).crctind)];
                    [IntExp(ExpIndex).deltadayCTlearning,IntExp(ExpIndex).deltanightCTlearning,IntExp(ExpIndex).AllDaysCTlearning,IntExp(ExpIndex).FFmorningfitCTlearning,IntExp(ExpIndex).FFeveningfitCTlearning]=daynight...
                        (IntExp(ExpIndex).FFdataCTlearning,IntExp(ExpIndex).fvalsCTlearning);    

                    IntExp(ExpIndex).FFdataFBbaseline=mean(experiment(1).contours(IntExp(ExpIndex).window,experiment(1).crfbind));
                    IntExp(ExpIndex).fvalsFBbaseline=experiment(1).fv(experiment(1).crfbind);
                    [IntExp(ExpIndex).deltadayFBbaseline,IntExp(ExpIndex).deltanightFBbaseline,IntExp(ExpIndex).AllDaysFBbaseline,IntExp(ExpIndex).FFmorningfitFBbaseline,IntExp(ExpIndex).FFeveningfitFBbaseline]=daynight...
                        (IntExp(ExpIndex).FFdataFBbaseline,IntExp(ExpIndex).fvalsFBbaseline);

                    IntExp(ExpIndex).FFdataCTbaseline=mean(experiment(1).contours(IntExp(ExpIndex).window,experiment(1).crctind));
                    IntExp(ExpIndex).fvalsCTbaseline=experiment(1).fv(experiment(1).crctind);
                    [IntExp(ExpIndex).deltadayCTbaseline,IntExp(ExpIndex).deltanightCTbaseline,IntExp(ExpIndex).AllDaysCTbaseline,IntExp(ExpIndex).FFmorningfitCTbaseline,IntExp(ExpIndex).FFeveningfitCTbaseline]=daynight...
                        (IntExp(ExpIndex).FFdataCTbaseline,IntExp(ExpIndex).fvalsCTbaseline);    
                % What is the reversion relative to baseline effect???
                            FBB=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysFBbaseline)
                                FBB=[FBB IntExp(ExpIndex).AllDaysFBbaseline(i).FFdata];
                            end
                            CTB=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysCTbaseline)
                                CTB=[CTB IntExp(ExpIndex).AllDaysCTbaseline(i).FFdata];
                            end
                            FBbaseFF=median(FBB);
                            CTbaseFF=median(CTB);
                            FBLearning=[];
                            CTLearning=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysFBlearning)
                                FBLearning(i)=median(IntExp(ExpIndex).AllDaysFBlearning(i).FFdata);
                                CTLearning(i)=median(IntExp(ExpIndex).AllDaysCTlearning(i).FFdata);
                            end
                            FBLearning-median(FBB)
                            CTLearning-median(CTB)
                            1-(FBLearning-median(FBB))./(CTLearning-median(CTB)) % Reversion %age (0.0 is no reversion, 1.0 is complete reversion)
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                         figure;hold on;
                        % baseline
                            i=1;
                            wide=5;
                            window=[300:350];
                                    a=[21:length(experiment(i).fv)];
                                    plot(runningaverage(timing4(experiment(i).fv(experiment(i).crfbind(experiment(i).crfbind>20))),wide),runningaverage(mean(experiment(i).contours(window,experiment(i).crfbind(experiment(i).crfbind>20))),wide),'r.','Markersize',15)
                                    plot(runningaverage(timing4(experiment(i).fv(experiment(i).crctind(experiment(i).crctind>20))),wide*9),runningaverage(mean(experiment(i).contours(window,experiment(i).crctind(experiment(i).crctind>20))),wide*9),'k.','Markersize',15)
                        % wn on  
                        for i=[2:3]% 7 14]
                                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                                    timedata=b;
                                    FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                                    lastpoints=find(diff(timedata)>10);
                                    startpoints=[1 lastpoints+1];
                                    endpoints=[lastpoints length(timedata)];
                                    for j=1:length(endpoints)
                                        if endpoints(j)>startpoints(j)+wide-1
                                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'r.','Markersize',15)
                                        end
                                    end
                                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                                    timedata=b;
                                    FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                                    lastpoints=find(diff(timedata)>10);
                                    startpoints=[1 lastpoints+1];
                                    endpoints=[lastpoints length(timedata)];
                                    for j=1:length(endpoints)
                                        if endpoints(j)>startpoints(j)+wide-1
                                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide*9),runningaverage(FFdata(startpoints(j):endpoints(j)),wide*9),'b.','Markersize',15)
                                        end

                                    end
                        end
                 %%%%%%%%%%%%%%%%
                 for i=1:length(IntExp(ExpIndex).FFmorningfitFBlearning)
                     t1=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(1));
                     t2=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(end));
                     plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBlearning(i);IntExp(ExpIndex).FFeveningfitFBlearning(i)],'k-','Linewidth',4)
                 end
                 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bk91w60
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    load /cardinal5/bk91w60/Experiment1019.mat
                   ExpIndex=3;
                    IntExp(ExpIndex).bird='bk91w60';    
                    IntExp(ExpIndex).window=[500:550];

                    IntExp(ExpIndex).FFdataFBlearning=[mean(experiment(4).contours(IntExp(ExpIndex).window,experiment(4).crfbind)) mean(experiment(5).contours(IntExp(ExpIndex).window,experiment(5).crfbind)) mean(experiment(6).contours(IntExp(ExpIndex).window,experiment(6).crfbind)) mean(experiment(7).contours(IntExp(ExpIndex).window,experiment(7).crfbind))];
                    IntExp(ExpIndex).fvalsFBlearning=[experiment(4).fv(experiment(4).crfbind) experiment(5).fv(experiment(5).crfbind) experiment(6).fv(experiment(6).crfbind) experiment(7).fv(experiment(7).crfbind)];
                    [IntExp(ExpIndex).deltadayFBlearning,IntExp(ExpIndex).deltanightFBlearning,IntExp(ExpIndex).AllDaysFBlearning,IntExp(ExpIndex).FFmorningfitFBlearning,IntExp(ExpIndex).FFeveningfitFBlearning]=daynight...
                        (IntExp(ExpIndex).FFdataFBlearning,IntExp(ExpIndex).fvalsFBlearning);


                    IntExp(ExpIndex).FFdataCTlearning=[mean(experiment(4).contours(IntExp(ExpIndex).window,experiment(4).crctind)) mean(experiment(5).contours(IntExp(ExpIndex).window,experiment(5).crctind)) mean(experiment(6).contours(IntExp(ExpIndex).window,experiment(6).crctind)) mean(experiment(7).contours(IntExp(ExpIndex).window,experiment(7).crctind))];
                    IntExp(ExpIndex).fvalsCTlearning=[experiment(4).fv(experiment(4).crctind) experiment(5).fv(experiment(5).crctind) experiment(6).fv(experiment(6).crctind) experiment(7).fv(experiment(7).crctind)];
                    [IntExp(ExpIndex).deltadayCTlearning,IntExp(ExpIndex).deltanightCTlearning,IntExp(ExpIndex).AllDaysCTlearning,IntExp(ExpIndex).FFmorningfitCTlearning,IntExp(ExpIndex).FFeveningfitCTlearning]=daynight...
                        (IntExp(ExpIndex).FFdataCTlearning,IntExp(ExpIndex).fvalsCTlearning);    


                    IntExp(ExpIndex).FFdataFBbaseline=[mean(experiment(2).contours(IntExp(ExpIndex).window,experiment(2).crfbind)) mean(experiment(3).contours(IntExp(ExpIndex).window,experiment(3).crfbind))];
                    IntExp(ExpIndex).fvalsFBbaseline=[experiment(2).fv(experiment(2).crfbind) experiment(3).fv(experiment(3).crfbind)];
                    [IntExp(ExpIndex).deltadayFBbaseline,IntExp(ExpIndex).deltanightFBbaseline,IntExp(ExpIndex).AllDaysFBbaseline,IntExp(ExpIndex).FFmorningfitFBbaseline,IntExp(ExpIndex).FFeveningfitFBbaseline]=daynight...
                        (IntExp(ExpIndex).FFdataFBbaseline,IntExp(ExpIndex).fvalsFBbaseline);


                    IntExp(ExpIndex).FFdataCTbaseline=[mean(experiment(2).contours(IntExp(ExpIndex).window,experiment(2).crctind)) mean(experiment(3).contours(IntExp(ExpIndex).window,experiment(3).crctind))];
                    IntExp(ExpIndex).fvalsCTbaseline=[experiment(2).fv(experiment(2).crctind) experiment(3).fv(experiment(3).crctind)];
                    [IntExp(ExpIndex).deltadayCTbaseline,IntExp(ExpIndex).deltanightCTbaseline,IntExp(ExpIndex).AllDaysCTbaseline,IntExp(ExpIndex).FFmorningfitCTbaseline,IntExp(ExpIndex).FFeveningfitCTbaseline]=daynight...
                        (IntExp(ExpIndex).FFdataCTbaseline,IntExp(ExpIndex).fvalsCTbaseline);    

                  % What is the reversion relative to baseline effect???
                            FBB=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysFBbaseline)
                                FBB=[FBB IntExp(ExpIndex).AllDaysFBbaseline(i).FFdata];
                            end
                            CTB=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysCTbaseline)
                                CTB=[CTB IntExp(ExpIndex).AllDaysCTbaseline(i).FFdata];
                            end
                            FBbaseFF=median(FBB);
                            CTbaseFF=median(CTB);
                            FBLearning=[];
                            CTLearning=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysFBlearning)
                                FBLearning(i)=median(IntExp(ExpIndex).AllDaysFBlearning(i).FFdata);
                                CTLearning(i)=median(IntExp(ExpIndex).AllDaysCTlearning(i).FFdata);
                            end
                            FBLearning-median(FBB)
                            CTLearning-median(CTB)
                            1-(FBLearning-median(FBB))./(CTLearning-median(CTB)) % Reversion %age (0.0 is no reversion, 1.0 is complete reversion)
  
 
            figure;hold on;
             wideCatches=[10 10 10 20 10 10 50 20];
             wideStims=[10 10 5 5 5 5 20 10];
            for i=3:7%length(experiment)
                        wideCatch=wideCatches(i);
                        wideStim=wideStims(i);
                        window=[500:550];
                                            [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                                timedata=b;
                                FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                                lastpoints=find(diff(timedata)>10);
                                startpoints=[1 lastpoints+1];
                                endpoints=[lastpoints length(timedata)];
                                for j=1:length(endpoints)
                                    if endpoints(j)>startpoints(j)+wideCatch-1
                                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wideCatch),runningaverage(FFdata(startpoints(j):endpoints(j)),wideCatch),'*','Color','k')
                                    end

                                end
                                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                                timedata=b;
                                FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                                lastpoints=find(diff(timedata)>10);
                                startpoints=[1 lastpoints+1];
                                endpoints=[lastpoints length(timedata)];
                                for j=1:length(endpoints)
                                    if endpoints(j)>startpoints(j)+wideStim-1
                                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wideStim),runningaverage(FFdata(startpoints(j):endpoints(j)),wideStim),'*','Color','r')
                                    end
                                end

            end                     
         %%%%%%%%%%%%%%%%
                 for i=1:length(IntExp(ExpIndex).FFmorningfitFBlearning)
                     t1=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(1));
                     t2=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(end));
                     plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBlearning(i);IntExp(ExpIndex).FFeveningfitFBlearning(i)],'k-','Linewidth',4)
                 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bk80w28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    load /bulbul1/LMAN_microstimulation//bk80w28/experimentB.mat
                   ExpIndex=4;
                    IntExp(ExpIndex).bird='bk80w28';    
                    IntExp(ExpIndex).window=[400:440];
                    
                    IntExp(ExpIndex).FFdataFBlearning=[];
                    IntExp(ExpIndex).fvalsFBlearning=[];                    
                    for i=1:9
                        IntExp(ExpIndex).FFdataFBlearning=[IntExp(ExpIndex).FFdataFBlearning mean(experiment(i+36).contours(IntExp(ExpIndex).window,experiment(i+36).crfbind))];
                        IntExp(ExpIndex).fvalsFBlearning=[IntExp(ExpIndex).fvalsFBlearning experiment(i+36).fv(experiment(i+36).crfbind)];
                    end
                    [IntExp(ExpIndex).deltadayFBlearning,IntExp(ExpIndex).deltanightFBlearning,IntExp(ExpIndex).AllDaysFBlearning,IntExp(ExpIndex).FFmorningfitFBlearning,IntExp(ExpIndex).FFeveningfitFBlearning]=daynight...
                        (IntExp(ExpIndex).FFdataFBlearning,IntExp(ExpIndex).fvalsFBlearning);
                    IntExp(ExpIndex).FFdataCTlearning=[];
                    IntExp(ExpIndex).fvalsCTlearning=[];
                    for i=1:9
                        IntExp(ExpIndex).FFdataCTlearning=[IntExp(ExpIndex).FFdataCTlearning mean(experiment(i+36).contours(IntExp(ExpIndex).window,experiment(i+36).crctind))];
                        IntExp(ExpIndex).fvalsCTlearning=[IntExp(ExpIndex).fvalsCTlearning experiment(i+36).fv(experiment(i+36).crctind)];
                    end
                    [IntExp(ExpIndex).deltadayCTlearning,IntExp(ExpIndex).deltanightCTlearning,IntExp(ExpIndex).AllDaysCTlearning,IntExp(ExpIndex).FFmorningfitCTlearning,IntExp(ExpIndex).FFeveningfitCTlearning]=daynight...
                        (IntExp(ExpIndex).FFdataCTlearning,IntExp(ExpIndex).fvalsCTlearning);    

                    IntExp(ExpIndex).FFdataFBbaseline=[mean(experiment(36).contours(IntExp(ExpIndex).window,experiment(36).crfbind))];
                    IntExp(ExpIndex).fvalsFBbaseline=[experiment(36).fv(experiment(36).crfbind)];
                    [IntExp(ExpIndex).deltadayFBbaseline,IntExp(ExpIndex).deltanightFBbaseline,IntExp(ExpIndex).AllDaysFBbaseline,IntExp(ExpIndex).FFmorningfitFBbaseline,IntExp(ExpIndex).FFeveningfitFBbaseline]=daynight...
                        (IntExp(ExpIndex).FFdataFBbaseline,IntExp(ExpIndex).fvalsFBbaseline);

                    IntExp(ExpIndex).FFdataCTbaseline=[mean(experiment(36).contours(IntExp(ExpIndex).window,experiment(36).crctind))];
                    IntExp(ExpIndex).fvalsCTbaseline=[experiment(36).fv(experiment(36).crctind) ];
                    [IntExp(ExpIndex).deltadayCTbaseline,IntExp(ExpIndex).deltanightCTbaseline,IntExp(ExpIndex).AllDaysCTbaseline,IntExp(ExpIndex).FFmorningfitCTbaseline,IntExp(ExpIndex).FFeveningfitCTbaseline]=daynight...
                        (IntExp(ExpIndex).FFdataCTbaseline,IntExp(ExpIndex).fvalsCTbaseline);    

                  % What is the reversion relative to baseline effect???
                            FBB=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysFBbaseline)
                                FBB=[FBB IntExp(ExpIndex).AllDaysFBbaseline(i).FFdata];
                            end
                            CTB=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysCTbaseline)
                                CTB=[CTB IntExp(ExpIndex).AllDaysCTbaseline(i).FFdata];
                            end
                            FBbaseFF=median(FBB);
                            CTbaseFF=median(CTB);
                            FBLearning=[];
                            CTLearning=[];
                            for i=1:length(IntExp(ExpIndex).AllDaysFBlearning)
                                FBLearning(i)=median(IntExp(ExpIndex).AllDaysFBlearning(i).FFdata);
                                CTLearning(i)=median(IntExp(ExpIndex).AllDaysCTlearning(i).FFdata);
                            end
                            FBLearning-median(FBB)
                            CTLearning-median(CTB)
                            1-(FBLearning-median(FBB))./(CTLearning-median(CTB)) % Reversion %age (0.0 is no reversion, 1.0 is complete reversion)
  
 
            figure;hold on;
             wideCatches=[10 10 10 20 10 10 50 20];
             wideStims=[10 10 5 5 5 5 20 10];
            for i=36:45%length(experiment)
                        wideCatch=0%wideCatches(i);
                        wideStim=0%wideStims(i);
                        window=[400:440];
                                            [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                                timedata=b;
                                FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                                lastpoints=find(diff(timedata)>10);
                                startpoints=[1 lastpoints+1];
                                endpoints=[lastpoints length(timedata)];
                                for j=1:length(endpoints)
                                    if endpoints(j)>startpoints(j)+wideCatch-1
                                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wideCatch),runningaverage(FFdata(startpoints(j):endpoints(j)),wideCatch),'o','markersize',5,'Color','k')
                                    end

                                end
                                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                                timedata=b;
                                FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                                lastpoints=find(diff(timedata)>10);
                                startpoints=[1 lastpoints+1];
                                endpoints=[lastpoints length(timedata)];
                                for j=1:length(endpoints)
                                    if endpoints(j)>startpoints(j)+wideStim-1
                                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wideStim),runningaverage(FFdata(startpoints(j):endpoints(j)),wideStim),'o','Color','r')
                                    end
                                end

            end                     
         %%%%%%%%%%%%%%%%
                 for i=1:length(IntExp(ExpIndex).FFmorningfitFBlearning)
                     t1=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(1));
                     t2=timing4(IntExp(ExpIndex).AllDaysFBlearning(i).fv(end));
                     plot([t1;t2],[IntExp(ExpIndex).FFmorningfitFBlearning(i);IntExp(ExpIndex).FFeveningfitFBlearning(i)],'k-','Linewidth',4)
                 end                 
                            
                            
                            
                            
                            
    
% <<Bird 1>> % bk48w74
        load /cardinal3/StimBirds/bk48w74/experiments.mat
        % bk48w74.m --> notes
        % bk48w74stim.m --> for plotting, etc.
        % Day/Night Intermittent
                % July 30th - August 4th - experiment(2) -- context is experiment(1:5)
                    % Learning --- Very good data, clear reversion and gradual consolidation
                                % only problem = 10% stim rate - undersampled
                    % Recovery --- Not ideal because consolidation never reached asymptote
                                % like like recovery in DMP happens overnight!
        % 100% stim during acquisition - acausal window
                % August 9th - experiment(6) -- context is experiment(5:8)
                % Looks like robust learning occurs
                % Note that experiment 7 stims throughout premotor window and shows
                     % reversion (although not original intent)
        % 100% stim during acquisition - acausal window
                % August 12th - experiment(10) -- context is experiment(8:12)*** caveat with exp 11 explained below
                % This is the same note and same direction as experiment(6)
                % Note that experiment(9) and experiment(10) are the same data - experiment 10 has more data though
                % Looks like robust learning occurs
              % AFTER THIS (Aug 13th) the computer crashed, oscillation, stim, fail
                % HOWEVER experiment 11 appears to show a bit of reversion (dubious)
                    % which could be used if necessary as evidence for stim working during
                    % experiment(10)
                % experiment(12) is nice for showing eventual recovery to baseline

                
% <<Bird 2>> % o9o4
        % o9o4.m
        % o9o4stim.m
            % Stopped singing when WN turned on (true before and after implantation)
            % Chronic stim
            % Problem ---> he occassionally/often pulled the lead out on one side of
                % his head in the morning ---> this was with the old leads/implants
            % 9.01-9.14 - dicking around with parameters, trying to get him to sing
        % Chronic stim experiment
            load Data917chronic.mat %????%
            % September 15-18 
            % no effect but pulled out probes intermittently
        % Chronic stim experiment
            load Data927.mat %????%
            % September 23-27
            % no effect


% <<Bird 3>> % bk91w60
        load /cardinal5/bk91w60/Experiment1019.mat
        % bk91w60.m --> notes
        % bk91w60stimanaly.m --> for plotting, etc.
        % Day/Night Intermittent
                % October 13th - October 18th - experiment(4:7) -- context is experiment(3:8)
                    % Learning --- Not much learning, some reversion but not that much, 
                                % rapid and then plateauing consolidation?
                                % only problem = 10% stim rate - undersampled
                    % Recovery --- Not ideal 
                  % In summary, the learning isn't sufficient magnitude (60Hz)
                  % and I think the stimulators are not sufficiently effective.
                  % However, it could be interpreted as a partial disruption.
             % NOTE that baseline stim effect appears to be upwards and if
             % you subtract this the results look better.
         % STIM CRAPPINESS CONFIRMED BY OSCILLOSCOPE   

         
% <<Bird 4>> % r87g80
    load /cardinal2/StimBirds/r87g80/ExperimentB_all.mat % second of two stack notes
     % COVERT LEARNING EXPERIMENT 1 with good positive controls --- upward=experiment(36:48) - November 11-19
        % Positive control - November 11-12th - experiment(37) --- context is experiment(36:38)
            % good learning - around 40Hz - threshold set at median+20Hz
        % 100% stim during acquisition - premotor window- November 14-15th - experiment(40) --- context is experiment(38:41)
            % no learning - threshold set at median+20Hz
        % 100% stim w/o WN - chronic stim - premotor window - November 16-17th - experiment(42) --- context is experiment(41:43)
            % no change
        % 100% stim during acquisition - acausal window - November 18-19th - experiment (46) --- context is experiment(44-48)
            % good learning - around 30-40Hz - threshold set at median+20Hz
            % note that in experiments 44-47 stim was in acausal window, but in experiment 48 stim was during the premotor window to cause reversion
            figure;hold on;xlim([360 700])
            plot(median(experiment(47).contours(:,experiment(47).crctind(1:50))')-median(experiment(44).contours(:,experiment(44).crctind(end-50:end))'),'b')
            plot(median(experiment(41).contours(:,experiment(41).crctind(1:50))')-median(experiment(38).contours(:,experiment(38).crctind(end-50:end))'),'k')
            plot(median(experiment(38).contours(:,experiment(38).crctind(1:50))')-median(experiment(36).contours(:,experiment(36).crctind(end-50:end))'),'r')
     % COVERT LEARNING EXPERIMENT 2 with good positive controls --- downward
        % Positive control - November 21-22 - experiment(49) --- context is experiment(48:50)
            % good learning - around 40Hz - threshold set at median-20Hz
      %%% Exps 51-61: Stimulator appeared to crap out and not block expression
      %%% Exp 62: Changed stim configs (November 29th)
        % 100% stim during acquisition - premotor window- December 1-2 - experiment(66) --- context is experiment(64-67)
            % minimal learning (maybe 15Hz) - threshold set at median-20Hz
        % 100% stim during acqusition - acausal window - December 3-4 - experiment(71) --- context is experiment(68:72)
       figure;hold on;xlim([360 700])    
       plot(median(experiment(50).contours(:,experiment(50).crctind(1:50))')-(median(experiment(48).contours(:,experiment(48).crctind(end-50:end))')),'b')                 
       plot(median(experiment(67).contours(:,experiment(67).crctind(1:50))')-(median(experiment(64).contours(:,experiment(64).crctind(end-50:end))')),'k')                 
       plot(median(experiment(72).contours(:,experiment(72).crctind(1:50))')-(median(experiment(68).contours(:,experiment(68).crctind(end-50:end))')),'r')     
     % COVERT LEARNING EXPERIMENTS with crappy positive controls --- experiment(18:35) - October 31-November 5
        % 100% stim during acqusition - premotor window -DOWN- October 31 - experiment(20) - context is experiment(18:21)
        % 100% stim during acqusition - premotor window -UP- November 1 - experiment(23) - context is experiment(21:24)          
        % 100% stim during acqusition - premotor window -DOWN-- November 2 - experiment(27) - context is experiment(25:28)
        % 100% stim during acqusition - premotor window -UP- November 3 - experiment(30) - context is experiment(28:31)
            figure;hold on;xlim([360 700])
            plot(-1*(median(experiment(21).contours(:,experiment(21).crctind(1:50))')-median(experiment(18).contours(:,experiment(18).crctind(end-50:end))')),'b')
            plot(median(experiment(24).contours(:,experiment(24).crctind)')-median(experiment(21).contours(:,experiment(21).crctind(end-50:end))'),'r')
            plot(-1*(median(experiment(28).contours(:,experiment(28).crctind(1:50))')-median(experiment(25).contours(:,experiment(25).crctind(end-50:end))')),'b')
            plot(median(experiment(31).contours(:,experiment(31).crctind(1:50))')-median(experiment(28).contours(:,experiment(28).crctind(end-50:end))'),'r') 
        % 100% stim w/o WN - chronic stim - premotor window - November 4 - experiment(32) - context is experiment(31:33)
        % 100% stim w/o WN - chronic stim - premotor window - November 4 - experiment(34) - context is experiment(33:35) 
            plot(median(experiment(33).contours(:,experiment(33).crctind(1:50))')-median(experiment(31).contours(:,experiment(31).crctind(end-50:end))'),'k')
            plot(median(experiment(35).contours(:,experiment(35).crctind(1:50))')-median(experiment(33).contours(:,experiment(33).crctind(end-50:end))'),'k')     
     % OTHER DATA - experiment(2:17) - October 19-27
        % This data may be useful at some point but probably not.
        % experiment(2-6) - 100% premotor stim+wn for 2-3 days leads to slight
            % increase in FF in same direction as stim offset and learning
            % direction
        % experiment(8-11) - 100% premotor stim w/o wn for 2-3 days
            % no effect, but stimulator crapped out
        % experiment(12-17) - stim was crapped out most of this period
    load /cardinal2/StimBirds/r87g80/experimentC.mat % first high stack after two low stack notes    
     % COVERT LEARNING
        % Positive control -DOWN-
            % Very slow
        % 100% stim during acquisition -UP- experiment(15) - context is experiment([12:17])
            % Definite increase in FF (~40Hz) but why? chronic stim or learning?
            % Why is positive control learning so slow?
            % Caveat of other high stack note having different mean FF? ---- learning to copy it?
        
        
% <<Bird 5>> % r30g38
    % November 2010 - high stack note
        % Tried stuff with high stack note but slow learning confounded things
    % December 2010 - low stack note
        % COVERT LEARNING
            % Positive control % experiment 3-4 (context is experiment 2-5)
            % 100% stim during acquisition % experiment 8 (context is experiment 6-9
                % stim effect started dying - failure to block learning
        % DAY/NIGHT INTERMITTENT % December 8-16
            % experiment 9 is pre; experiment 10:11 is wn on (pooled together which includes 12&13); experiment 14 is wn off
            % note that during experiment 14 I accidentally shorted one of
            %    the stim sides - so ignore recovery data - see r30g38.m notes
            % Note somewhat of a concern about incomplete inactivation/reversion
    % January 2010
        % problems - inability to completely block expression with 100% stim

% <<Bird 6>> % bk80w28
    % December 2010
    % COVERT LEARNING
        % Positive control % experiment 8 (context is experiment 7-10)
        % 100% stim during acquistion - premotor window % experiment 14 (context is experiment 12-15)
            % failure to block learning despite offset (but no s.d. reduction)
    % January 2010
    % problems with reduction in stim rate with 100% stim
    % attempt reversion experiment with 10-20% stim

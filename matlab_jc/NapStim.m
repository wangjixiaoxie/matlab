

%%% 06.09.11
% Add Tim's data
% Only look at first 4 days
% consider outliers

% Step 1: evaluate code in NapStimpre1.m
clearvars -except ExpJC           
% OR...
load /bulbul1/LMAN_microstimulation/NapStim/processedExps.mat
% 1. syll time, 2. isFB?, 3. FF, 4. song index, 5. syll index

                % Only consider first four days
                for i=1:length(ExpJC)
                    if max(ExpJC(i).performance(:,1))-min(ExpJC(i).performance(:,1))>48
                        nights=find(diff(ExpJC(i).performance(:,1))>8);
                    else
                        nights=find(diff(ExpJC(i).performance(:,1))>0.3);
                    end
                    if length(nights)>3
                        ExpJC(i).performance=ExpJC(i).performance(1:nights(4),:);
                    end
                end
                %
                % Remove outliers using the 1.0 IQR rule (very conservative - 1.5 is
                % typical)
                % ExpJCoutlier=ExpJC;
                for i=1:length(ExpJC)
                    runmdn=runningmedian(ExpJC(i).performance(:,3),100);
                    runprchigh=runningprctile(ExpJC(i).performance(:,3),75,100);
                    runprclow=runningprctile(ExpJC(i).performance(:,3),25,100);
                    runprchigh2=[ones(1,50)*runprchigh(1) runprchigh ones(1,50)*runprchigh(end)];
                    runprclow2=[ones(1,50)*runprclow(1) runprclow ones(1,50)*runprclow(end)];   
                    runiqrhigh=runprchigh2+1*(runprchigh2-runprclow2);
                    runiqrlow=runprclow2-1*(runprchigh2-runprclow2);
                    % column 6 - isOutlier?
                    for j=1:length(ExpJC(i).performance(:,3))
                        ExpJC(i).performance(j,6)=0;
                        if ExpJC(i).performance(j,3)>runiqrhigh(j) | ExpJC(i).performance(j,3)<runiqrlow(j)
                            ExpJC(i).performance(j,6)=1;
                        end
                    end
                   ExpJC(i).performance=ExpJC(i).performance(~ExpJC(i).performance(:,6),:);
                end
                %



                %
                %
                for i=1:length(ExpJC)
                    performance=ExpJC(i).performance;
                    indFB=find(performance(:,2));
                    indCT=find(~performance(:,2));
                    clear stimtrials
                    stimtrials=performance(indFB,:);
                    catchtrials=performance(indCT,:);
                    [e,f]=unique(stimtrials(:,4));
                    count=0;
                    for j=1:length(f)
                        d=find(stimtrials(:,4)==e(j));
                        for k=1:length(d)
                            count=count+1;
                            stimtrials(count,6)=k;
                        end
                    end

                    % Correlation between catch & stim
                    count=0;
                    clear mnFFStim mnFFCatch
                    for j=1:max(performance(:,4)) % for each song
                        clear indS indC
                        indS=find(performance(:,4)==j & performance(:,2)==1);
                        indC=find(performance(:,4)==j & performance(:,2)==0);
                        if ~isempty(indS) & ~isempty(indC)
                            count=count+1;
                            mnFFStim(count)=mean(performance(indS,3));
                            mnFFCatch(count)=mean(performance(indC,3));
                        end
                    end
                    var1=diff(mnFFStim);
                    var2=diff(mnFFCatch);
                    a=corrcoef(var1,var2); % 0.0174
                    stimVScatch(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));


                    %
                    % Created a vector in syllable order with times (1), FB/CT (2), FF (3),
                    % song number (4), and count of stimmed syllable in that specific song (5)


                    % NAP/GAP-DEPENDENT CONSOLIDATION
                    % 1 - Correlation between total time between consecutive stimulations
                    % and consolidation
                    var1=diff(stimtrials(:,1));
                    var2=diff(stimtrials(:,3));
                    a=corrcoef(var1,var2); % 0.0169 n.s.
                    TimeBT_all(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                    % ignoring nights
                    ind1=find(var1<8);
                    var1=var1(ind1);
                    var2=var2(ind1);
                    a=corrcoef(var1,var2); % 0.0247 n.s.
                    TimeBT_nonight(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                    % collapsing stims in one song onto a single point
                    clear stimsongtime stimsongFF stimsongind
                    [e,f]=unique(stimtrials(:,4));
                    for j=1:length(f)
                        % mean time of song
                        stimsongtime(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),1));
                        % mean FF of song
                        stimsongFF(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),3));
                        % number of song
                        stimsongind(j)=e(j);
                    end
                    var1=diff(stimsongtime);
                    var2=diff(stimsongFF);
                    a=corrcoef(var1,var2); % 0.0174
                    TimeBT_songcollapse(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                    % 2 - Correlation between gap and consolidation in cases with stim in
                    % consecutive songs.
                    indconsecutive=find(diff(stimsongind)==1);
                    var1=var1(indconsecutive);
                    var2=var2(indconsecutive);
                    a=corrcoef(var1,var2); % 0.0174
                    TimeBT_songcollapse_consec(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2))); % 0.0232



                    % PERFORMANCE-DEPENDENT CONSOLIDATION
                    % 4 Difference between consolidation of multiple stim notes in a
                    % single song
                    var1=stimtrials(:,6); % position in song (as # of stimmed note)
                    var2=stimtrials(:,3); % FF
                    a=corrcoef(var1,var2); % -0.1205 --- forgetting in a single song???
                    positioncorr(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                    % alternatively, to show it isn't just a bias (e.g.) of songs with few
                    % stims happening late in learning
                    % Compare 2nd vs. 1st
                    a1=find(var1==2);
                    MN_one_back(i)=mean(var2(a1)-var2(a1-1));
                    % Compare 3rd vs 2nd vs 1st
                    a1=find(var1==3);
                    mean(var2(a1)-var2(a1-1)); % -7.90
                    mean(var2(a1)-var2(a1-2)); % -18.88
                    % BUT...Could this also be happening at baseline????

                    % 5 Correlation between number of renditions sung between two stimmed syllables and consolidation
                    var1=diff(stimtrials(:,5)); % syllables since most recent
                    var2=diff(stimtrials(:,3)); % consolidation since most recent
                    a=corrcoef(var1,var2); %
                    performancecorr(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                    % Collapse across songs
                    clear stimsongperf stimsongFF
                    [e,f]=unique(stimtrials(:,4));
                    for j=1:length(f)
                        stimsongperf(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),5));
                        stimsongFF(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),3));
                    end
                    var1=diff(stimsongperf);
                    var2=diff(stimsongFF);
                    a=corrcoef(var1,var2);
                    performancecorr_songs(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));

                    % 6 Correlation between number of renditions sung between first

                end


                % All values should be positive
                % All values are t-statistics
                % Naptime
                goodind=[2:7];
                TimeBT_all % 1. Correlation between total time between consecutive stimulations and consolidation in that period
                TimeBT_nonight % 2. Same as 1, but ignoring overnight periods
                TimeBT_songcollapse % 3. Same as 1, but collapses across songs with multiple stim trials
                TimeBT_songcollapse_consec % 4. Same as 4, but only considers consecutive songs 
                % that have stim, to avoid caveat that most of the time was consumed by singing (as opposed to naps)

                % Practice makes perfect
                positioncorr % 5. In some songs, multiple syllables are stimmed.  Does position in song correlate with consolidation? 
                MN_one_back % 6. Is consolidation in the 2nd stimmed syllable of song greater than in the first? (similar to 5)
                performancecorr % 7. Correlation between total number of syllables sung between two stimmed syllables and consolidation?
                performancecorr_songs % significant negative correlation % 8. Same as 7, but collapses across songs (to avoid caveat of #5-6)

                % Okay wow, none of these things are good predictors.
                % Another possible predictor is that the amount of learning that occurs in
                % catch notes
                % Do this by song to avoid confounds:
                stimVScatch; % 11. 

                % CONCLUSIONS
                    % When there are several non-stimmed songs between stimmed songs, you get
                    % less consolidation




                % Look at baseline - I've only done this for the first 4 experiments (mine)

                            % Get 9-10
                            clear positioncorrPRE MN_one_backPRE
                            for i=1:length(IntEXp) % For each experiment
                                % Create a vector in syllable order with times (1), FB/CT (2), FF (3), and song
                                % number (4)
                                clear songtimeCT sylltimeCT songtimeFB sylltimeFB syllFFCT syllFFFB
                                for j=1:length(IntExp(i).fvalsCTbaseline)
                                    syllFFCT(j)=IntExp(i).FFdataCTbaseline(j);
                                    songtimeCT(j)=timing4(IntExp(i).fvalsCTbaseline(j));
                                    sylltimeCT(j)=timing4(IntExp(i).fvalsCTbaseline(j))+(1/1000)*(1/3600)*IntExp(i).fvalsCTbaseline(j).ons(IntExp(i).fvalsCTbaseline(j).ind);
                                end
                                for j=1:length(IntExp(i).fvalsFBbaseline)
                                    syllFFFB(j)=IntExp(i).FFdataFBbaseline(j);
                                    songtimeFB(j)=timing4(IntExp(i).fvalsFBbaseline(j));
                                    sylltimeFB(j)=timing4(IntExp(i).fvalsFBbaseline(j))+(1/1000)*(1/3600)*IntExp(i).fvalsFBbaseline(j).ons(IntExp(i).fvalsFBbaseline(j).ind);
                                end
                                isempty(find(diff(sylltimeFB<0))) % should be true if in order
                                isempty(find(diff(sylltimeCT<0))) % should be true if in order
                                isempty(find(diff(sylltimeFB==0))) % should be true if in order
                                isempty(find(diff(sylltimeCT==0))) % should be true if in order
                                songsboth=[songtimeCT songtimeFB];
                                timesboth=[sylltimeCT sylltimeFB];
                                FFboth=[syllFFCT syllFFFB];
                                [a,b]=sort(timesboth);
                                FFsort=FFboth(b);
                                songsort=songsboth(b);
                                [c,d]=unique(songsort);
                                clear performance
                                for j=1:length(timesboth)
                                    if b(j)>length(sylltimeCT)
                                        performance(j,2)=1; % FB
                                    else
                                        performance(j,2)=0; % CT
                                    end
                                    performance(j,1)=a(j);
                                    performance(j,3)=IntExp(i).dir*FFsort(j); % normalized for baseline direction
                                    performance(j,4)=find(songsort(j)==c);
                                    performance(j,5)=j;
                                end
                                indFB=find(performance(:,2));
                                indCT=find(~performance(:,2));
                                clear stimtrials
                                stimtrials=performance(indFB,:);
                                catchtrials=performance(indCT,:);
                                [e,f]=unique(stimtrials(:,4));
                                count=0;
                                for j=1:length(f)
                                    d=find(stimtrials(:,4)==e(j));
                                    for k=1:length(d)
                                        count=count+1;
                                        stimtrials(count,6)=k;
                                    end
                                end

                                % PERFORMANCE-DEPENDENT CONSOLIDATION
                                % 4 Difference between consolidation of multiple stim notes in a
                                % single song
                                var1=stimtrials(:,6); % position in song (as # of stimmed note)
                                var2=stimtrials(:,3); % FF
                                a=corrcoef(var1,var2); % -0.1205 --- forgetting in a single song???
                                positioncorrPRE(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                                % alternatively, to show it isn't just a bias (e.g.) of songs with few
                                % stims happening late in baseline
                                % Compare 2nd vs. 1st
                                a1=find(var1==2);
                                MN_one_backPRE(i)=mean(var2(a1)-var2(a1-1));
                                % Compare 3rd vs 2nd vs 1st
                                a1=find(var1==3);
                                mean(var2(a1)-var2(a1-1)); % -7.90
                                mean(var2(a1)-var2(a1-2)); % -18.88
                                % BUT...Could this also be happening at baseline????
                            end


                positioncorrPRE % 9. Baseline --- In some songs, multiple syllables are stimmed.  Does position in song correlate with consolidation? 
                MN_one_backPRE % 10. Baseline --- Is consolidation in the 2nd stimmed syllable of song greater than in the first? (similar to 5)



    


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%


            % 06.02.11

            load /cardinal7/LMAN_microstimulation/IntermittentAnaly060211.mat

                        clear TimeBT_all TimeBT_nonight TimeBT_songcollapse TimeBT_songcollapse_consec positioncorr MN_one_back performancecorr performancecorr_songs
                        for i=1:length(IntExp) % For each experiment
                            % Create a vector in syllable order with times (1), FB/CT (2), FF (3), and song
                            % number (4)
                                    clear songtimeCT sylltimeCT songtimeFB sylltimeFB syllFFCT syllFFFB
                                    for j=1:length(IntExp(i).fvalsCTlearning)
                                        syllFFCT(j)=IntExp(i).FFdataCTlearning(j);
                                        songtimeCT(j)=timing4(IntExp(i).fvalsCTlearning(j));
                                        sylltimeCT(j)=timing4(IntExp(i).fvalsCTlearning(j))+(1/1000)*(1/3600)*IntExp(i).fvalsCTlearning(j).ons(IntExp(i).fvalsCTlearning(j).ind);        
                                    end    
                                    for j=1:length(IntExp(i).fvalsFBlearning)
                                        syllFFFB(j)=IntExp(i).FFdataFBlearning(j);       
                                        songtimeFB(j)=timing4(IntExp(i).fvalsFBlearning(j));
                                        sylltimeFB(j)=timing4(IntExp(i).fvalsFBlearning(j))+(1/1000)*(1/3600)*IntExp(i).fvalsFBlearning(j).ons(IntExp(i).fvalsFBlearning(j).ind);
                                    end
                                    isempty(find(diff(sylltimeFB<0))) % should be true if in order
                                    isempty(find(diff(sylltimeCT<0))) % should be true if in order
                                    isempty(find(diff(sylltimeFB==0))) % should be true if in order
                                    isempty(find(diff(sylltimeCT==0))) % should be true if in order
                                    songsboth=[songtimeCT songtimeFB];
                                    timesboth=[sylltimeCT sylltimeFB];
                                    FFboth=[syllFFCT syllFFFB];
                                    [a,b]=sort(timesboth);
                                    FFsort=FFboth(b);
                                    songsort=songsboth(b);
                                    [c,d]=unique(songsort);
                                    clear performance
                                    for j=1:length(timesboth)
                                        if b(j)>length(sylltimeCT)
                                            performance(j,2)=1; % FB
                                        else
                                            performance(j,2)=0; % CT
                                        end
                                        performance(j,1)=a(j);
                                        performance(j,3)=IntExp(i).dir*FFsort(j); % normalized for learning direction
                                        performance(j,4)=find(songsort(j)==c);
                                        performance(j,5)=j;
                                    end
                                    indFB=find(performance(:,2));
                                    indCT=find(~performance(:,2));
                                    clear stimtrials
                                    stimtrials=performance(indFB,:);
                                    catchtrials=performance(indCT,:);
                                    [e,f]=unique(stimtrials(:,4));
                                    count=0;
                                    for j=1:length(f)
                                        d=find(stimtrials(:,4)==e(j));
                                        for k=1:length(d)
                                            count=count+1;
                                            stimtrials(count,6)=k;
                                        end
                                    end

                                % Correlation between catch & stim
                                count=0;
                                clear mnFFStim mnFFCatch
                                for j=1:max(performance(:,4)) % for each song
                                    clear indS indC
                                    indS=find(performance(:,4)==j & performance(:,2)==1);
                                    indC=find(performance(:,4)==j & performance(:,2)==0);
                                    if ~isempty(indS) & ~isempty(indC)
                                        count=count+1;
                                        mnFFStim(count)=mean(performance(indS,3));
                                        mnFFCatch(count)=mean(performance(indC,3));
                                    end
                                end
                                var1=diff(mnFFStim);
                                var2=diff(mnFFCatch);
                                a=corrcoef(var1,var2); % 0.0174
                                stimVScatch(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));


                                %
                        % Created a vector in syllable order with times (1), FB/CT (2), FF (3),
                        % song number (4), and count of stimmed syllable in that specific song (5)


                        % NAP/GAP-DEPENDENT CONSOLIDATION            
                             % 1 - Correlation between total time between consecutive stimulations
                             % and consolidation
                                        var1=diff(stimtrials(:,1));
                                        var2=diff(stimtrials(:,3));
                                        a=corrcoef(var1,var2); % 0.0169 n.s.
                                        TimeBT_all(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                                    % ignoring nights
                                        ind1=find(var1<8);
                                        var1=var1(ind1);
                                        var2=var2(ind1);
                                        a=corrcoef(var1,var2); % 0.0247 n.s.  
                                        TimeBT_nonight(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                                    % collapsing stims in one song onto a single point
                                        clear stimsongtime stimsongFF stimsongind
                                        [e,f]=unique(stimtrials(:,4));
                                        for j=1:length(f)
                                            % mean time of song
                                            stimsongtime(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),1));
                                            % mean FF of song
                                            stimsongFF(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),3));
                                            % number of song
                                            stimsongind(j)=e(j);
                                        end
                                        var1=diff(stimsongtime);
                                        var2=diff(stimsongFF);
                                        a=corrcoef(var1,var2); % 0.0174
                                        TimeBT_songcollapse(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                             % 2 - Correlation between gap and consolidation in cases with stim in
                             % consecutive songs.
                                        indconsecutive=find(diff(stimsongind)==1);
                                        var1=var1(indconsecutive);
                                        var2=var2(indconsecutive);
                                        a=corrcoef(var1,var2); % 0.0174
                                        TimeBT_songcollapse_consec(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2))); % 0.0232



                        % PERFORMANCE-DEPENDENT CONSOLIDATION
                             % 4 Difference between consolidation of multiple stim notes in a
                             % single song
                              var1=stimtrials(:,6); % position in song (as # of stimmed note)
                              var2=stimtrials(:,3); % FF
                              a=corrcoef(var1,var2); % -0.1205 --- forgetting in a single song???
                              positioncorr(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                                % alternatively, to show it isn't just a bias (e.g.) of songs with few
                                % stims happening late in learning
                                    % Compare 2nd vs. 1st
                                        a1=find(var1==2);
                                        MN_one_back(i)=mean(var2(a1)-var2(a1-1));
                                    % Compare 3rd vs 2nd vs 1st
                                        a1=find(var1==3);
                                        mean(var2(a1)-var2(a1-1)); % -7.90
                                        mean(var2(a1)-var2(a1-2)); % -18.88
                                     % BUT...Could this also be happening at baseline????

                             % 5 Correlation between number of renditions sung between two stimmed syllables and consolidation
                              var1=diff(stimtrials(:,5)); % syllables since most recent
                              var2=diff(stimtrials(:,3)); % consolidation since most recent
                              a=corrcoef(var1,var2); % 
                              performancecorr(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                                % Collapse across songs
                                        clear stimsongperf stimsongFF
                                        [e,f]=unique(stimtrials(:,4));
                                        for j=1:length(f)
                                            stimsongperf(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),5));
                                            stimsongFF(j)=mean(stimtrials(find(stimtrials(:,4)==e(j)),3));
                                        end
                                        var1=diff(stimsongperf);
                                        var2=diff(stimsongFF);
                                        a=corrcoef(var1,var2); 
                                        performancecorr_songs(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));

                              % 6 Correlation between number of renditions sung between first 

                        end


                        % Get 9-10
                        clear positioncorrPRE MN_one_backPRE
                        for i=1:length(IntExp) % For each experiment
                            % Create a vector in syllable order with times (1), FB/CT (2), FF (3), and song
                            % number (4)
                                    clear songtimeCT sylltimeCT songtimeFB sylltimeFB syllFFCT syllFFFB
                                    for j=1:length(IntExp(i).fvalsCTbaseline)
                                        syllFFCT(j)=IntExp(i).FFdataCTbaseline(j);
                                        songtimeCT(j)=timing4(IntExp(i).fvalsCTbaseline(j));
                                        sylltimeCT(j)=timing4(IntExp(i).fvalsCTbaseline(j))+(1/1000)*(1/3600)*IntExp(i).fvalsCTbaseline(j).ons(IntExp(i).fvalsCTbaseline(j).ind);        
                                    end    
                                    for j=1:length(IntExp(i).fvalsFBbaseline)
                                        syllFFFB(j)=IntExp(i).FFdataFBbaseline(j);       
                                        songtimeFB(j)=timing4(IntExp(i).fvalsFBbaseline(j));
                                        sylltimeFB(j)=timing4(IntExp(i).fvalsFBbaseline(j))+(1/1000)*(1/3600)*IntExp(i).fvalsFBbaseline(j).ons(IntExp(i).fvalsFBbaseline(j).ind);
                                    end
                                    isempty(find(diff(sylltimeFB<0))) % should be true if in order
                                    isempty(find(diff(sylltimeCT<0))) % should be true if in order
                                    isempty(find(diff(sylltimeFB==0))) % should be true if in order
                                    isempty(find(diff(sylltimeCT==0))) % should be true if in order
                                    songsboth=[songtimeCT songtimeFB];
                                    timesboth=[sylltimeCT sylltimeFB];
                                    FFboth=[syllFFCT syllFFFB];
                                    [a,b]=sort(timesboth);
                                    FFsort=FFboth(b);
                                    songsort=songsboth(b);
                                    [c,d]=unique(songsort);
                                    clear performance
                                    for j=1:length(timesboth)
                                        if b(j)>length(sylltimeCT)
                                            performance(j,2)=1; % FB
                                        else
                                            performance(j,2)=0; % CT
                                        end
                                        performance(j,1)=a(j);
                                        performance(j,3)=IntExp(i).dir*FFsort(j); % normalized for baseline direction
                                        performance(j,4)=find(songsort(j)==c);
                                        performance(j,5)=j;
                                    end
                                    indFB=find(performance(:,2));
                                    indCT=find(~performance(:,2));
                                    clear stimtrials
                                    stimtrials=performance(indFB,:);
                                    catchtrials=performance(indCT,:);
                                    [e,f]=unique(stimtrials(:,4));
                                    count=0;
                                    for j=1:length(f)
                                        d=find(stimtrials(:,4)==e(j));
                                        for k=1:length(d)
                                            count=count+1;
                                            stimtrials(count,6)=k;
                                        end
                                    end

                        % PERFORMANCE-DEPENDENT CONSOLIDATION
                             % 4 Difference between consolidation of multiple stim notes in a
                             % single song
                              var1=stimtrials(:,6); % position in song (as # of stimmed note)
                              var2=stimtrials(:,3); % FF
                              a=corrcoef(var1,var2); % -0.1205 --- forgetting in a single song???
                              positioncorrPRE(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                                % alternatively, to show it isn't just a bias (e.g.) of songs with few
                                % stims happening late in baseline
                                    % Compare 2nd vs. 1st
                                        a1=find(var1==2);
                                        MN_one_backPRE(i)=mean(var2(a1)-var2(a1-1));
                                    % Compare 3rd vs 2nd vs 1st
                                        a1=find(var1==3);
                                        mean(var2(a1)-var2(a1-1)); % -7.90
                                        mean(var2(a1)-var2(a1-2)); % -18.88
                                     % BUT...Could this also be happening at baseline????
                        end
            % All values should be positive
            % All values are t-statistics
            % Naptime
            TimeBT_all % 1. Correlation between total time between consecutive stimulations and consolidation in that period
            TimeBT_nonight % 2. Same as 1, but ignoring overnight periods
            TimeBT_songcollapse % 3. Same as 1, but collapses across songs with multiple stim trials
            TimeBT_songcollapse_consec % 4. Same as 4, but only considers consecutive songs 
            % that have stim, to avoid caveat that most of the time was consumed by singing (as opposed to naps)

            % Practice makes perfect
            positioncorr % 5. In some songs, multiple syllables are stimmed.  Does position in song correlate with consolidation? 
            MN_one_back % 6. Is consolidation in the 2nd stimmed syllable of song greater than in the first? (similar to 5)
            performancecorr % 7. Correlation between total number of syllables sung between two stimmed syllables and consolidation?
            performancecorr_songs % 8. Same as 7, but collapses across songs (to avoid caveat of #5-6)
            positioncorrPRE % 9. Baseline --- In some songs, multiple syllables are stimmed.  Does position in song correlate with consolidation? 
            MN_one_backPRE % 10. Baseline --- Is consolidation in the 2nd stimmed syllable of song greater than in the first? (similar to 5)


            % Okay wow, none of these things are good predictors.
            % Another possible predictor is that the amount of learning that occurs in
            % catch notes
            % Do this by song to avoid confounds:
            stimVScatch % 11. 

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
        % 06.01.11
        %
        load /cardinal7/LMAN_microstimulation/IntermittentAnaly041211b.mat
        %
        % CATCH
        tbase=timing4(IntExp(4).fvalsCTbaseline);
        tlearn=timing4(IntExp(4).fvalsCTlearning);
        tbaseCT=tbase-min(tlearn);
        tlearnCT=tlearn-min(tlearn);
        FFbaseCT=IntExp(4).FFdataCTbaseline-mean(IntExp(4).FFdataCTbaseline);
        FFlearnCT=IntExp(4).FFdataCTlearning-mean(IntExp(4).FFdataCTbaseline);
        figure;hold on;
        plot(tbaseCT,FFbaseCT,'k.')
        plot(tlearnCT,FFlearnCT,'r.')

        corrcoef(tlearnCT,FFlearnCT) % -0.6047
        corrcoef([1:1:length(FFlearnCT)],FFlearnCT) % -0.6076

        % FB
        indthis=4;
        tbase=timing4(IntExp(4).fvalsFBbaseline);
        tlearn=timing4(IntExp(4).fvalsFBlearning);
        tbaseFB=tbase-min(tlearn);
        tlearnFB=tlearn-min(tlearn);
        FFbaseFB=IntExp(4).FFdataFBbaseline-mean(IntExp(4).FFdataFBbaseline);
        FFlearnFB=IntExp(4).FFdataFBlearning-mean(IntExp(4).FFdataFBbaseline);
        figure;hold on;
        plot(tbaseFB,FFbaseFB,'k.')
        plot(tlearnFB,FFlearnFB,'b.')

        %

        %
        % For each FB sung
        clear etime indACT learnA learnC gappvs
        for i=2:length(FFlearnFB)
            nextone=min(find(tlearnCT>tlearnFB(i)));   
            pvsone=max(find(tlearnCT<tlearnFB(i)));       
            if isempty(nextone);return;end    
            % What is the elapsed time?
            etime(i)=tlearnFB(i);
            % How much time since most recent song?
            gappvs(i)=etime(i)-tlearnCT(pvsone);
            % How many times has the syllable been sung;
            indACT(i)=nextone;
            % What is the amount of actual learning?
            learnA(i)=mean(FFlearnCT(nextone-1:nextone));
            % How much consolidation?
            learnC(i)=FFlearnFB(i);
        end
        % correlations with learnC
            % etime R=-0.5701
            % learnC R=0.4980
            % indACT R=-0.5737

            [b,bint,r] = regress(learnC',[etime;learnA;indACT]'); 
            % significant effects of learnA and indACT but not etime

        % Is there a significant effect of time indep of syllables sung?
            [b1,bint1,r1] = regress(learnC',indACT'); 
            [b2,bint2,r2] = regress(r1,etime'); % no effect
        % Is there a significant effect of syllables sung indep of time?
            [b3,bint3,r3] = regress(learnC',etime'); 
            [b4,bint4,r4] = regress(r3,indACT'); % no effect

        % Is there any correlation between amount of time between performances and
        % consolidation at that time?
        var1=diff(etime);
        var2=learnC(2:end);
        corrcoef(var1,var2) % NO - r=-0.0334
        % What about if you ignore nights
        indnap=find(var1<10);
        corrcoef(var1(indnap),var2(indnap)) % Still NO - r=-0.0365

        % Is there any correlation between number of performances between stimmed performances and
        % consolidation at that time?
            var1=diff(etime);
            var2=learnC(2:end);
            corrcoef(var1,var2) % NO - r=-0.0334
            % What about if you ignore nights
            indnap=find(var1<10);
            corrcoef(var1(indnap),var2(indnap)) % Still NO - r=-0.0365

        % To do: look at true gaps - gaps between pvs song and current song (not just between stimmed notes)
            corrcoef(gappvs,learnC) % NO - r=-0.0443
            % What about if you ignore nights?
            indnap=find(var1<10);
            corrcoef(gappvs(indnap),learnC(indnap)) % still no - r=-0.0445
                        % But this is actually close to significance - directional
                        % P=0.0622

        % For each song...
        [t1]=unique(tlearnFB);
        clear mnFFsong naptime
        for i=2:length(t1)
            % time of that song
                t1(i);
            % stimmed syllables in that song
                indstim=find(tlearnFB==t1(i));
            % average FF of those syllables
                mnFFsong(i)=mean(FFlearnFB(indstim));
            % time since most recent song
                pvsone=max(find(tlearnCT<tlearnFB(i)));       
                pvstime=tlearnCT(pvsone);
                naptime(i)=t1(i)-pvstime;
        end
        % What is the correlation between amnt of consolidation b/t songs and amnt
        % of time b/t those songs?
        corrcoef(diff(mnFFsong),naptime(2:end)) % R=0.0012 - no relationship - 
            % I checked this to make sure the indices were correct!





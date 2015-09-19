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


                    %
                    % Created a vector in syllable order with times (1), FB/CT (2), FF (3),
                    % song number (4), and count of stimmed syllable in that specific song (5)


                    % NAP/GAP-DEPENDENT CONSOLIDATION
                    % 1 - Correlation between total time between consecutive stimulations
                    % and consolidation
                    var1=diff(catchtrials(:,1));
                    var2=diff(catchtrials(:,3));
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
                    [e,f]=unique(catchtrials(:,4));
                    for j=1:length(f)
                        % mean time of song
                        stimsongtime(j)=mean(catchtrials(find(catchtrials(:,4)==e(j)),1));
                        % mean FF of song
                        stimsongFF(j)=mean(catchtrials(find(catchtrials(:,4)==e(j)),3));
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
                    var1=catchtrials(:,6); %+ position in song (as # of stimmed note)
                    var2=catchtrials(:,3); % FF
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
                    var1=diff(catchtrials(:,5)); % syllables since most recent
                    var2=diff(catchtrials(:,3)); % consolidation since most recent
                    a=corrcoef(var1,var2); %
                    performancecorr(i)=a(2)*sqrt((length(var1)-2)/(1-a(2)*a(2)));
                    % Collapse across songs
                    clear stimsongperf stimsongFF
                    [e,f]=unique(catchtrials(:,4));
                    for j=1:length(f)
                        stimsongperf(j)=mean(catchtrials(find(catchtrials(:,4)==e(j)),5));
                        stimsongFF(j)=mean(catchtrials(find(catchtrials(:,4)==e(j)),3));
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
                TimeBT_all; % 1. Correlation between total time between consecutive stimulations and consolidation in that period
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

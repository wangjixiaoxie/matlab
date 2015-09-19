load /bulbul1/LMAN_microstimulation/IntermittentAnaly060211.mat
load /stimsumoutJC.mat % TW - from dropbox
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
                        performance(:,1)=performance(:,1)/24;
                        ExpJC(i).performance=performance;
            end
            
            % Convert TW into my coordinates
            ExpTW=stimsumout;
            ljc=length(ExpJC);
            for i=1:length(ExpTW)
                clear performance
                for j=1:size(ExpTW(i).CTlearnvls,1)
                    ExpTW(i).CTlearnvls(j,5)=0;
                end
                for j=1:size(ExpTW(i).STlearnvls,1)
                    ExpTW(i).STlearnvls(j,5)=1;
                end
                allvls=[ExpTW(i).CTlearnvls;ExpTW(i).STlearnvls];
                [a,b]=sort(allvls(:,3)); % a are the times, b are the indices
                allsort=allvls(b,:);
                performance(:,1)=allsort(:,3); % syll time
                performance(:,2)=allsort(:,5); % Is FB?
                performance(:,3)=allsort(:,2)*0.5; % FF
                [c,d]=unique(allvls(:,1));
                for j=1:length(allsort)
                    performance(j,4)=(find(c==allsort(j,1)));
                    performance(j,5)=j;
                end
                ExpJC(ljc+i).performance=performance;
        
            end
    ExpJC(7).performance(:,3)=mean(stimsumout(3).CTbasvls(:,2)*0.5)-1*(ExpJC(7).performance(:,3)-mean(stimsumout(3).CTbasvls(:,2)*0.5));
            % 

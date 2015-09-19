% TW muscimol
% 11.08.11

%%%%%%%%%%%
% BIRD: PK20R49
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=1; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs


beginwn=avls.wn(1).on(1);

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
     plot(avls.adjvls{1}{1,i}(:,1)-beginwn,avls.adjvls{1}{1,i}(:,2),'.')
     if ~isempty([find(i==avls.muon)])
         i
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
     end
end
figure;hold on;
for i=1:length(avls.adjvls{2}) % for each run
     plot(avls.adjvls{2}{1,i}(:,1)-beginwn,avls.adjvls{2}{1,i}(:,2),'.')
     if ~isempty([find(i==avls.muon)])
         i
            plot(avls.adjvls{2}{1,i}(find(avls.adjvls{2}{1,i}(:,3)),1)-beginwn,avls.adjvls{2}{1,i}(find(avls.adjvls{2}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{2}{1,i}(find(1-avls.adjvls{2}{1,i}(:,3)),1)-beginwn,avls.adjvls{2}{1,i}(find(1-avls.adjvls{2}{1,i}(:,3)),2),'r.') % non-targeted
     end
end

% do these by hand
beginwn=avls.wn(1).on(1);
endwn=avls.wnrev(1).on(1);

% acsf - positive control - Hz/day
        tvacsf=[]; % normalized
        fvacsf=[]; % normalized
        allFFs=[]; % absolute
        alltimes=[]; % absolute
        for i=1:length(avls.adjvls{1}) % for each run
            if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if ACSF (i.e. not MUSC)
                if min(avls.adjvls{1}{1,i}(:,1))>beginwn & max(avls.adjvls{1}{1,i}(:,1))<endwn % if it is a wn run (i.e. not baseline or reversion)
                    i
                    if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                        tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                        FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                        tvacsf=[tvacsf tplus'];
                        fvacsf=[fvacsf FFplus'];
                        allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                        alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                    else
                        lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                        daybegins=[1 lightsout+1];
                        dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                        for j=1:length(lightsout) % for each day
                            tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                            FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                            tvacsf=[tvacsf tplus'];
                            fvacsf=[fvacsf FFplus'];
                            allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                            alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                        end
                    end
                end
            end
        end
        % check to make sure it looks reasonable
            [vl,ind1]=sort(tvacsf);
            figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
        % measure the slope in Hz/day
        [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
        slope=p(1); 
        mnFF=mean(allFFs);
        
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 

        alltimesort=alltimes(ind2);
        allFFsort=allFFs(ind2);
        
% muscimol
        tvmu=[];
        fvmu=[];
        %mnFFs=[];
        count=0;
        clear tvalspre tvalspost fvalspre fvalspost
        for i=1:length(avls.adjvls{1}) % for each run
            if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if MUSC
                onset=min(avls.adjvls{1}{1,i}(:,1));
                offset=max(avls.adjvls{1}{1,i}(:,1));
                if onset>beginwn & offset<endwn % if it is a wn run (i.e. not baseline or reversion)
                    count=count+1;
                    i
                    reltimespre=alltimesort-onset;
                    reltimespost=alltimesort-offset;
                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;  
                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                    fvalspre{count}=allFFsort(find(reltimespre<0));
                    fvalspost{count}=allFFsort(find(reltimespost>0));
                end
            end
        end
timegap=3/14;
% save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

% post-processing
clear Tpre FFpre Npre Tpost FFpost Npost
for i=1:length(tvalspre)
    Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
    FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
    CVpre(i)=std(fvalspre{i}(find(tvalspre{i}*24>-8)))/mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
    Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
    Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
    FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
    CVpost(i)=std(fvalspost{i}(find(tvalspost{i}*24>-8)))/mean(fvalspost{i}(find(tvalspost{i}*24>-8)));    
    Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
end
                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap
                CVpost
                CVpre

    
    
    
    
%%%%%%%%%%%
% BIRD: PU34
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=2; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(2);
beginwn=avls.wn(1).on(1);
load pathvals1-analdata.mat
figure;hold on;

for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        plot(avls.adjvls{1}{1,i}(:,1)-beginwn,avls.adjvls{1}{1,i}(:,2),'.')
        if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 6];
daymusc=[5 6];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
%                         figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=3.5/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost CVpre CVpost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                CVpre(i)=std(fvalspre{i}(find(tvalspre{i}*24>-8)))/mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                CVpost(i)=std(fvalspost{i}(find(tvalspost{i}*24>-8)))/mean(fvalspost{i}(find(tvalspost{i}*24>-8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
            timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
            NULL=timediffAC*slope*(14/24) % null model
            ACTUAL=FFpost-FFpre
            POSCTL=slope*timegap
            CVpost
            CVpre

% Experiment 2
dayacsf=[12 19];
daymusc=[16 17];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
            timegap=3.8666/14
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost CVpre CVpost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                CVpre(i)=std(fvalspre{i}(find(tvalspre{i}*24>-8)))/mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                CVpost(i)=std(fvalspost{i}(find(tvalspost{i}*24>-8)))/mean(fvalspost{i}(find(tvalspost{i}*24>-8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
            timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
            NULL=timediffAC*slope*(14/24) % null model
            ACTUAL=FFpost-FFpre
            POSCTL=slope*timegap
            CVpost
            CVpre
    
                
 %%%%%%%%%%%
% BIRD: BK20BK45
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=3; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(3);
beginwn=avls.wn(1).on(1);
load pathvals2-analdata.mat

figure;hold on;

for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        plot(avls.adjvls{1}{1,i}(:,1)-beginwn,avls.adjvls{1}{1,i}(:,2),'.')
        if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end
figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 10];
daymusc=[4 5];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=4/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost CVpre CVpost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=median(fvalspre{i}(find(tvalspre{i}*24>-8)));
                CVpre(i)=jcstd(fvalspre{i}(find(tvalspre{i}*24>-8)))/mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=median(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                CVpost(i)=jcstd(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)))/mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
            timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
            NULL=timediffAC*slope*(14/24) % null model
            ACTUAL=FFpost-FFpre
            POSCTL=slope*timegap
            CVpre
            CVpost
                
 %%%%%%%%%%%
% BIRD: BK20BK45
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=4; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(4);
beginwn=avls.wn(1).on(1);
load pathvals1-analdata.mat

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 10];
daymusc=[8 9];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=4/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap
 %%%%%%%%%%%
% BIRD: BK61W42
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=5; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(i);
beginwn=avls.wn(1).on(1);
load pathvals1-analdata.mat

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 3];
daymusc=[2 3];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=3.75/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap

                
                
  %%%%%%%%%%%
% BIRD: BK28W6
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=6; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(i);
beginwn=avls.wn(1).on(1);
load pathvals1-analdata.mat

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 3];
daymusc=[2 3];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=4.7/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap
               
  %%%%%%%%%%%
% BIRD: PK32BK28
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=7; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(i);
beginwn=avls.wn(1).on(1);
load pathvals3-analdata.mat

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 8];
daymusc=[6 8];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=3/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap
                
                
  %%%%%%%%%%%
% BIRD: BK61W42
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=8; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(i);
beginwn=avls.wn(1).on(1);
load pathvals3-analdata.mat

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

  %%%%%%%%%%%
% BIRD: PU57W52
% all musc runs have some escapes
%%%%%%%%%%%
clear all
edit birdstructlist % evaluate it - creates bs
edit shiftanal12
i=10; % for each index in bs
[sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, i)
% load the mat file specified in bs
sumbs=sumbs(i);
beginwn=avls.wn(1).on(1);
load pathvals1-analdata.mat

figure;hold on;
for i=1:length(avls.adjvls{1}) % for each run
    if ~isempty(avls.adjvls{1}{1,i}(:,1))
        if ~isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{2})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'b.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'g.') % non-targeted
        end
        if isempty([find(i==sumbs.shiftruns{1}) find(i==sumbs.shiftruns{1})])
            plot(avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(avls.adjvls{1}{1,i}(:,3)),2),'k.') % targeted
            plot(avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),1)-beginwn,avls.adjvls{1}{1,i}(find(1-avls.adjvls{1}{1,i}(:,3)),2),'r.') % non-targeted
        end
    end
end

% Experiment 1
dayacsf=[0 8];
daymusc=[6 8];

            % acsf - positive control - Hz/day
                    tvacsf=[]; % normalized
                    fvacsf=[]; % normalized
                    allFFs=[]; % absolute
                    alltimes=[]; % absolute
                    for i=1:length(avls.adjvls{1}) % for each run
                        if isempty([find(i==sumbs.shiftruns{1})])    % if ACSF (i.e. not MUSC)
                            if min(avls.adjvls{1}{1,i}(:,1))>beginwn+dayacsf(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+dayacsf(2) % if it is a wn run (i.e. not baseline or reversion)
                                if max(diff(avls.adjvls{1}{1,i}(:,1)))<0.33  % if there are no nights
                                    tplus=avls.adjvls{1}{1,i}(:,1)-min(avls.adjvls{1}{1,i}(:,1));
                                    FFplus=avls.adjvls{1}{1,i}(:,2)-mean(avls.adjvls{1}{1,i}(:,2));
                                    tvacsf=[tvacsf tplus'];
                                    fvacsf=[fvacsf FFplus'];
                                    allFFs=[allFFs avls.adjvls{1}{1,i}(:,2)'];
                                    alltimes=[alltimes avls.adjvls{1}{1,i}(:,1)'];
                                else
                                    lightsout=find(diff(avls.adjvls{1}{1,i}(:,1))>0.33); % returns last index b/f night of sleep
                                    daybegins=[1 lightsout+1];
                                    dayends=[lightsout length(avls.adjvls{1}{1,i}(:,1))];
                                    for j=1:length(lightsout) % for each day
                                        tplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)-min(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1));
                                        FFplus=avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)-mean(avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2));
                                        tvacsf=[tvacsf tplus'];
                                        fvacsf=[fvacsf FFplus'];
                                        allFFs=[allFFs avls.adjvls{1}{1,i}(daybegins(j):dayends(j),2)'];
                                        alltimes=[alltimes avls.adjvls{1}{1,i}(daybegins(j):dayends(j),1)'];
                                    end
                                end
                            end
                        end
                    end
                    % check to make sure it looks reasonable
                        [vl,ind1]=sort(tvacsf);
                        figure;plot(runningaverage(tvacsf(ind1),20),runningaverage(fvacsf(ind1),20))
                    % measure the slope in Hz/day
                    [p]=polyfit(tvacsf(ind1),fvacsf(ind1),1);
                    slope=p(1); 
                    mnFF=mean(allFFs);

                    %
        % THIS REMOVES SLEEP FROM + CTL CALCULATION (BETTER)
        [vl,ind2]=sort(alltimes);
        gaps=find(diff(alltimes(ind2))>0.3);
        alltimesawake=alltimes(ind2);
        for i=1:length(gaps)
            alltimesawake(gaps(i)+1:end)=alltimesawake(gaps(i)+1:end)-10/24;
        end
        p=polyfit(alltimesawake-min(alltimesawake),allFFs(ind2),1);
        slope=p(1); 
                    alltimesort=alltimes(ind2);
                    allFFsort=allFFs(ind2);

            % muscimol
                    tvmu=[];
                    fvmu=[];
                    %mnFFs=[];
                    count=0;
                    clear tvalspre tvalspost fvalspre fvalspost
                    for i=1:length(avls.adjvls{1}) % for each run
                        if min(avls.adjvls{1}{1,i}(:,1))>beginwn+daymusc(1) & max(avls.adjvls{1}{1,i}(:,1))<beginwn+daymusc(2) % if it is the run I'm looking for
                            if ~isempty([find(i==sumbs.shiftruns{1})])    % if MUSC
                                onset=min(avls.adjvls{1}{1,i}(:,1));
                                offset=max(avls.adjvls{1}{1,i}(:,1));
                                    count=count+1;
                                    i
                                    reltimespre=alltimesort-onset;
                                    reltimespost=alltimesort-offset;
                                    tvalspre{count}=alltimesort(find(reltimespre<0))-onset;
                                    tvalspost{count}=alltimesort(find(reltimespost>0))-offset;
                                    fvalspre{count}=allFFsort(find(reltimespre<0));
                                    fvalspost{count}=allFFsort(find(reltimespost>0));
                                    avls.tmoff(i)
                                    avls.tmon(i)
                            end
                        end
                    end
timegap=3/14;
            % save tvalspre, tvalspost, fvalspre, fvalspost, slope, mnFF

            % post-processing
            clear Tpre FFpre Npre Tpost FFpost Npost
            for i=1:length(tvalspre)
                Tpre(i)=mean(tvalspre{i}(find(tvalspre{i}*24>-8)));
                FFpre(i)=mean(fvalspre{i}(find(tvalspre{i}*24>-8)));
                Npre(i)=length((tvalspre{i}(find(tvalspre{i}*24>-8))));
                Tpost(i)=mean(tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                FFpost(i)=mean(fvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8)));
                Npost(i)=length((tvalspost{i}(find(tvalspost{i}*24>2 & tvalspost{i}*24<8))));
            end
                timediffAC=Tpost-Tpre; % How much time elapses with ACSF present? Since we can't get an instantaneous measure of FF
                NULL=timediffAC*slope*(14/24) % null model
                ACTUAL=FFpost-FFpre
                POSCTL=slope*timegap

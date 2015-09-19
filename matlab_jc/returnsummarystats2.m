function Exp1=returnsummarystats(Exp)
% Called by RepeatAnalysisNovember.m
% For each data set, it takes each day of singing and calculates repeat
% length, transition probabilities, hit rate (if applicable), and error
% bars for each of these statistics.

for i=1:length(Exp)
    i
    clearvars -except Exp Exp1 i
    if Exp(i).learning % For learning experiments
        MINI=Exp(i).MIN;
        MINIless=MINI-1;
        targetedtrans=MINIless:1:20; % targeted transitions
        % SPLIT ALL repeat length data into days
        daycount=0;
        for j=1:length(Exp(i).time) % Pre, Dur, Post
            if Exp(i).JC==0 % Time data for Evren vs. JC
                unidays=find(diff(Exp(i).time{j})>0.3); % Evren time
            else
                unidays=find(diff(Exp(i).time{j})>8); % JC time
            end
            unistart=[1 unidays+1];
            uniend=[unidays length(Exp(i).time{j})];
            for ii=1:length(unistart) % For each day
                daycount=daycount+1;
                allrpdata{daycount}=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                group(daycount)=j;
            end
        end
        % For each day, figure out all parameters
        count=0;
        for day=1:length(allrpdata) % for each day
            rpdata=allrpdata{day};
            % Repeat length & hit rate
            replong(day)=mean(rpdata);
            hitrate(day)=sum(rpdata(find(rpdata>MINIless))-MINIless*ones(1,length(find(rpdata>MINIless))))/sum(rpdata);
            % Create surrogate data sets using resampling
            clear replongBS hitrateBS
            for bs=1:1000
                scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                scrambleRPL=rpdata(scrambleIND);
                replongBS(bs)=mean(scrambleRPL);
                hitrateBS(bs)=sum(scrambleRPL(find(scrambleRPL>MINIless))-MINIless*ones(1,length(find(scrambleRPL>MINIless))))/sum(scrambleRPL);
            end
            % The standard deviation of these surrogate data sets is equivalent
            % to the standard error of the mean for the statistic of
            % interest
            semreplong(day)=std(replongBS);
            semhitrate(day)=std(hitrateBS);
            % Transition probabilities
            for k=1:1:20 % for each transition
                ptrans(day,k)=length(find(rpdata>k))/length(find(rpdata>k-1));
                samplesizeptrans(day,k)=length(find(rpdata>k));
                clear ptransBS
                for bs=1:1000
                    scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                    scrambleRPL=rpdata(scrambleIND);
                    ptransBS(bs)=length(find(scrambleRPL>k))/length(find(scrambleRPL>k-1));
                end
                semptrans(day,k)=std(ptransBS);
            end
        end
    else
        targetedtrans=[]; % targeted transitions
        % SPLIT ALL repeat length data into days
        daycount=0;
        for j=1 % Pre only
            if Exp(i).JC==0 % Time data for Evren vs. JC
                nightlength=0.3333; % Evren
            else
                nightlength=8; % JC
            end
            unidays=find(diff(Exp(i).time{j})>nightlength); 
            unistart=[1 unidays+1];
            uniend=[unidays length(Exp(i).time{j})];
            for ii=1:length(unistart) % For each day
                daycount=1+ceil((min(Exp(i).time{j}(unistart(ii):uniend(ii)))-min(Exp(i).time{j}))/(3*nightlength));
                allrpdata{ii}=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                group(daycount)=1;
            end
        end
        % For each day, figure out all parameters
        count=0;
        for day=1:length(allrpdata) % for each day
            rpdata=allrpdata{day};
            % Repeat length & hit rate
            replong(day)=mean(rpdata);
            % Create surrogate data sets using resampling
            clear replongBS hitrateBS
            for bs=1:1000
                scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                scrambleRPL=rpdata(scrambleIND);
                replongBS(bs)=mean(scrambleRPL);
            end
            % The standard deviation of these surrogate data sets is equivalent
            % to the standard error of the mean for the statistic of
            % interest
            semreplong(day)=std(replongBS);
            % Transition probabilities
            for k=1:1:20 % for each transition
                ptrans(day,k)=length(find(rpdata>k))/length(find(rpdata>k-1));
                samplesizeptrans(day,k)=length(find(rpdata>k));
                clear ptransBS
                for bs=1:1000
                    scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                    scrambleRPL=rpdata(scrambleIND);
                    ptransBS(bs)=length(find(scrambleRPL>k))/length(find(scrambleRPL>k-1));
                end
                semptrans(day,k)=std(ptransBS);
            end
        end
        hitrate=[];
        semhitrate=[];
    end
    
    Exp1(i).bird=Exp(i).bird;
    Exp1(i).hitbelow=Exp(i).hitbelow;
    Exp1(i).postlesion=Exp(i).postlesion;
    Exp1(i).learning=Exp(i).learning;
    Exp1(i).baselineShort=Exp(i).baselineShort;
    Exp1(i).baselineLong=Exp(i).baselineLong;
    Exp1(i).predays=find(group==1); % Indices for pre days
    Exp1(i).wndays=find(group==2); % Indices for wn days
    Exp1(i).postdays=find(group==3); % Indices for postdays
    Exp1(i).replong=replong; % Repeat lengths: replong(day,length)
    Exp1(i).hitrate=hitrate; % Hit rate: hitrate(day,rate)
    Exp1(i).semreplong=semreplong;
    Exp1(i).semhitrate=semhitrate;
    Exp1(i).targetedtrans=targetedtrans; % Which ptrans were targeted?
    Exp1(i).ptrans=ptrans; % ptrans(day,trans#)
    Exp1(i).samplesizeptrans=samplesizeptrans;
    Exp1(i).semptrans=semptrans;
end

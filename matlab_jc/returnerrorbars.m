function Exp=returnerrorbars(Exp)
% June 20, 2010 - RepeatAnalysisNovember

% no lesion, learning to shift down
% For each experiment, calculate ptrans, repeat lengths, and hit rates
for i=ind1;
    MINI=Exp(i).MIN;
    MINIless=MINI-1;
    clear lastptrans replong hitrate day
    count=0;
    clear semptrans semreplong semhitrate
    for j=1:3 % Pre, Dur, Post
            if Exp(i).JC==0
                unidays=find(diff(Exp(i).time{j})>0.3); % Evren time
            else
                unidays=find(diff(Exp(i).time{j})>8); % JC time
            end
            unistart=[1 unidays+1];
            uniend=[unidays length(Exp(i).time{j})];
            if i==6 & j==3 % Missing data for this experiment (Evren)
                lastptrans(j).data(1)=0;
                day(j).data(1)=0;
                replong(j).data(1)=0;
                hitrate(j).data(1)=0;
                semptrans(j).data(1)=0;
                semreplong(j).data(1)=0;
                semhitrate(j).data(1)=0;
            else

            for ii=1:length(unistart) % For each day
                rpdata=Exp(i).rplength{j}(unistart(ii):uniend(ii));
                count=count+1;
                lastptrans(j).data(ii)=length(find(rpdata>=MINI))/length(find(rpdata>=MINIless));
                day(j).data(ii)=count;
                replong(j).data(ii)=mean(rpdata);
                hitrate(j).data(ii)=sum(rpdata(find(rpdata>MINIless))-MINIless*ones(1,length(find(rpdata>MINIless))))/sum(rpdata);
                % Create surrogate data sets using resampling
                clear ptransBS replongBS hitrateBS
                for bs=1:1000
                    scrambleIND=ceil(rand(1,length(rpdata))*length(rpdata));
                    scrambleRPL=rpdata(scrambleIND);
                    ptransBS(bs)=length(find(scrambleRPL>=MINI))/length(find(scrambleRPL>=MINIless));
                    replongBS(bs)=mean(scrambleRPL);
                    hitrateBS(bs)=sum(scrambleRPL(find(scrambleRPL>MINIless))-MINIless*ones(1,length(find(scrambleRPL>MINIless))))/sum(scrambleRPL);
                end
                % The standard deviation of these surrogate data sets is equivalent
                % to the standard error of the mean for the statistic of
                % interest
                semptrans(j).data(ii)=std(ptransBS);
                semreplong(j).data(ii)=std(replongBS);
                semhitrate(j).data(ii)=std(hitrateBS);
            end
        end
    end
    [probstay,numstay,semstay]=repeatsTOptrans2(Exp,i);
    % Resample these probabilities

    Exp(i).probstay=probstay;
    Exp(i).numstay=numstay;
    Exp(i).semstay=semstay;
    Exp(i).lastptrans=lastptrans;
    Exp(i).day=day;
    Exp(i).replong=replong;
    Exp(i).hitrate=hitrate;
    Exp(i).semptrans=semptrans;
    Exp(i).semreplong=semreplong;
    Exp(i).semhitrate=semhitrate;
end

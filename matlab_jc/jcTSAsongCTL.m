function [CTL]=jcTSAsongCTL(batch,pitch,startingpoint)

for ii=1:size(pitch,2)
    FFstats(ii)=median(pitch(startingpoint-150:startingpoint,ii));
end

[vals,trigs]=triglabel(batch,'a',1,1,0,0);
aa=zeros(length(trigs),8);
count=0;
for i=1:length(trigs)-4
    for j=1:length(trigs(i).ntmintmp)
        count=count+1;
        if FFstats(count)>prctile(FFstats,55)
            aa(i,j)=1;
        else
            aa(i,j)=2;
        end
    end
end

counttot=0;
countESC=0;
samesong=[];
nextsong=[];
ones=[];
twos=[];
threes=[];
fours=[];
nextplace=[];
for k=1:length(trigs)
    for kk=1:7
        if aa(k,kk)~=0
            counttot=counttot+1;
        end
        if aa(k,kk)==1 %if escape
            countESC=countESC+1;
            escapes(countESC)=FFstats(counttot);
            addon=0;
            for ll=kk+1:8 %the rest of the song
                addon=addon+1;
                if aa(k,ll)~=0
                    samesong=[samesong FFstats(counttot+addon)];
                else totsong=ll-1;
                    break
                end
            end
            addnext=0;
            if k~=length(trigs)
                for lll=1:7 %the next song
                    addnext=addnext+1;
                    if aa(k+1,lll)~=0
                        nextsong=[nextsong FFstats(counttot+addnext+totsong-kk)];
                    end
                end
            end
            if k~=length(trigs)
                if aa(k+1,kk)~=0 %Same position in the next song
                    nextplace=[nextplace FFstats(counttot+totsong)];
                end
            end
            ones=[ones FFstats(counttot+1)];
            twos=[twos FFstats(counttot+2)];
            threes=[threes FFstats(counttot+3)];
            fours=[fours FFstats(counttot+4)];
        end
    end
end
CTL.medsame=median(samesong);
CTL.medplace=median(nextplace);
CTL.mednext=median(nextsong);
CTL.medFF=median(FFstats);
CTL.medesc=median(escapes);
CTL.serr=std(nextsong)/sqrt(length(nextsong));
CTL.medones=median(ones);
CTL.medtwos=median(twos);
CTL.medthrees=median(threes);
CTL.medfours=median(fours);
plot(1,CTL.medFF,'*','Color','r');hold on;
plot(1.1,CTL.medsame,'*')
plot(1.4,CTL.medplace,'*')
plot(1.7,CTL.mednext,'*')
plot(2,CTL.medesc,'*','Color','g')
plot(2.1,CTL.medones,'*')
plot(2.2,CTL.medtwos,'*')
plot(2.3,CTL.medthrees,'*')
plot(2.4,CTL.medfours,'*')
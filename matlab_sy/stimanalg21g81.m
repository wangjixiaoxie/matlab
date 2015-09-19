%master script, load and combine various fvs, from postscreen./50microstim,
%and 10 microstim

pathvl{1}='/doyale2/g21g81/10microstim/'
catchvl{1}='batch.catch'
pathvl{2}='/doyale2/g21g81/50microstim/'
catchvl{2}='batch.keep.catch'
pathvl{3}='/doyale2/g21g81/postscreen/'
catchvl{3}='batch.keep.comb'
pathvl{4}='/doyale1/twarren/g21g81/stimoff/'
catchvl{4}='batch.keep.catch'


fvcomb=[];
makefv=1;
calc_prop_note=0;
%stimstart 1/30 12pm

if (makefv)
    for ii=1:length(pathvl)
        NT='a';PRENT='';PSTNT='';
        tbinshft=0.012;
        strcmd=['cd '  pathvl{ii}];
        eval(strcmd);
        fvnam='fva'
        bt=catchvl{ii}
        NFFT=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
        fbins=[2500,4500; 5000,9000];
        strcmd=[fvnam '=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,''obs0'');']
        eval(strcmd);
        strcmd=['save -append ' bt '.mat ' fvnam    ' tbinshft ' ' fbins']
        eval(strcmd)
    end
end


%2/2 1355 - switch to make template earlier.
%2/5 1635 - switch to make cutoff higher.

clear switchdts

bnds{1}='2007-01-22 07:00:00'
bnds{2}='2007-03-07 07:00:00'

%start, <3900
switchdt{1}='2007-01-30 12:00:00';
%trigger time earlier
switchdt{2}='2007-02-02 13:55:00';
%
switchdt{3}='2007-02-05 16:35:00';
%
switchdt{4}='2007-02-07 20:00:00'
%
switchdt{5}='2007-02-13 13:35:00'
%stim off
switchdt{6}='2007-02-15 19:00:00'
for ii=1:length(switchdt)
    switchdts(ii)=datenum(switchdt{ii},'yyyy-mm-dd HH:MM:SS');
end
for ii=1:length(bnds)
    bndsjs(ii)=datenum(bnds{ii},'yyyy-mm-dd HH:MM:SS');
end

fvcomb=[];
for ii=1:length(pathvl)
    strcmd=['load ' pathvl{ii} catchvl{ii} '.mat']
    eval(strcmd);
    fvtmp=fva;
    fvcomb=[fvcomb fvtmp];
end
vals=getvals(fvcomb,2,'TRIG');
vlot=ceil((vals(:,1)-bndsjs(1)));
switchdays=reduce_dates(switchdts,bndsjs(1));
daysbnds=[bndsjs(2)-bndsjs(1)]
figure
hold on;
days=[1:daysbnds]
unique(vlot);
calc_prop_note=0;
for ii=1:length(days)
    indtmp=find(vlot==days(ii))
    mnvl(ii)=mean(vals(indtmp,2))
    stdv(ii)=std(vals(indtmp,2))
    nnotes(ii)=length(indtmp)
    if (calc_prop_note)
        svfnval=fvcomb(indtmp(1))
        targntcnt=0;
        totntcnt=0;
        for jj=1:length(indtmp)
            currfnval=fv(indtmp(jj)).fn;
            if(currfnval==svfnval)
                targntcnt=ntcnt+1;
                tempcnt=fv(indtmp(jj)).lbl;
            elseif jj==length(indtmp)
                totntcnt=totntcnt+tmpcnt
            else
                svfnval=curfnval;
                totntcnt=totntcnt+tmpcnt;
            end
        end
propnote(ii)=targntcnt/totntcnt;
     end
end


figure

x=[switchdays switchdays];
y=[4200 4400];
plot(x,y,'r','Linewidth',3)

hold on;
errorbar(days, mnvl, stdv,'k+','Linewidth',3)





%plotting hit rates
%distributions of hits and misses

indhit=find(vals(:,3)==1)
indmiss=find(vals(:,3)==0)
daystm=unique(vlot(indhit));

for ii=1:length(daystm)
    indyht=find(vlot(indhit)==daystm(ii));
    indyms=find(vlot(indmiss)==daystm(ii));
    htrt(ii)=length(indyht)/(length(indyms)+length(indyht));
end
figure
plot(daystm,htrt,'-o');
hold on;
y=[.8 1]
plot(x,y,'r','Linewidth',3)

%below plots total number of hits on each day.
figure
x=[switchdays switchdays]
y=[1000 9000]
plot(x,y,'r', 'Linewidth',3)
hold on;
newdt=datevec(vals(:,1));
plot(newdt(:,3), vals(:,3),'o')


%analysis of syntax for b.
%probability of transition

bttl=[]
bsepttl=[]
conttl=[]

for ii=1:length(days)
    bcnt=0;
    bsepcnt=0;
    concnt=0;
    
    dyvl=days(ii);
    ind=find(vlot==dyvl)
    for jj=1:length(ind)
        indvl=ind(jj);
        if(fvcomb(indvl).trans(1)=='b')
            bcnt=bcnt+1;
        elseif (fvcomb(indvl).trans=='-b')
            bsepcnt=bsepcnt+1;
    
        else
            concnt=concnt+1;
        end
    end
    bttl(dyvl)=bcnt;
    bsepttl(dyvl)=bsepcnt;
    conttl(dyvl)=concnt;
end

for ii=1:dyvl
    
    sumvl(ii)=bttl(ii)+bsepttl(ii)+conttl(ii)
    if(sumvl(ii)==0)
       bpct(ii)=0
       bsepct(ii)=0
       conpct(ii)=0
    
    else
        bpct(ii)=bttl(ii)/sumvl(ii)
        bsepct(ii)=bsepttl(ii)/sumvl(ii)
        conpct(ii)=conttl(ii)/sumvl(ii)
    end
end

figure;
hold on;
plot(days, bpct, 'bo')
plot(days,bsepct,'ro')
plot(days,conpct,'ko')

attl=[]
endttl=[]
conttl=[]
bttl=[]

%analysis of syntax for a
for ii=1:length(days)
    acnt=0;
    endcnt=0;
    concnt=0;
    bcnt=0;
    dyvl=days(ii);
    
    ind=find(vlot==dyvl)
    for jj=1:length(ind)
        indvl=ind(jj)
        if(fvcomb(indvl).trans=='-a')
            acnt=acnt+1;
        elseif (fvcomb(indvl).trans=='nd')
            endcnt=endcnt+1;
        elseif (fvcomb(indvl).trans=='ii')
            endcnt=endcnt+1;
        elseif (fvcomb(indvl).trans(1)=='b'|fvcomb(indvl).trans(2)=='b')
            bcnt=bcnt+1;
        else
            concnt=concnt+1;
        end
    end
    attl(dyvl)=acnt;
    endttl(dyvl)=endcnt
    conttl(dyvl)=concnt
    bttl(dyvl)=bcnt
    
    for ii=1:dyvl
    
        sumvl(ii)=bttl(ii)+endttl(ii)+conttl(ii)+attl(ii)
        if(sumvl(ii)==0)
            apct(ii)=0
            endpct(ii)=0
            conpct(ii)=0
            bpct(ii)=0
        else
        bpct(ii)=bttl(ii)/sumvl(ii)
        endpct(ii)=bsepttl(ii)/sumvl(ii)
        conpct(ii)=conttl(ii)/sumvl(ii)
        apct(ii)=attl(ii)/sumvl(ii)
        end
    end
end

figure
hold on
plot(days,apct,'ko')
plot(days,bpct,'ro')
plot(days,conpct,'go')
plot(days,endpct,'bo')


    




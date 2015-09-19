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
calc_prop_note=1;

if (makefv)
    
    for ii=1:length(pathvl)
    
        NT='b';PRENT='';PSTNT='';
        tbinshft=0.012;
        strcmd=['cd '  pathvl{ii}];
        eval(strcmd);
        
        fvnam='fvb'
        bt=catchvl{ii}
        
        NFFT=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
        fbins=[2500,4500; 5000,9000];
       % edges=[6000:75:8000];
    
        strcmd=[fvnam '=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,''obs0'');']
        
        eval(strcmd);
        
        strcmd=['save ' bt '.mat ' fvnam    ' tbinshft ' ' fbins']
        eval(strcmd)
        
    end
    
    %vals=getvals(fv,2,'TRIG');

end

%stimstart 1/30 12pm

%2/2 1355 - switch to make template earlier.
%2/5 1635 - switch to make cutoff higher.

clear switchdt

bnds{1}='2007-01-25 07:00:00'
bnds{2}='2007-03-07 07:00:00'


switchdt{1}='2007-01-30 12:00:00';
switchdt{2}='2007-02-02 13:55:00';
switchdt{3}='2007-02-05 16:35:00';
switchdt{4}='2007-02-07 20:00:00'
switchdt{5}='2007-02-13 13:35:00'
switchdt{6}='2007-02-15 19:00:00'
calc_prop_note=0;

for ii=1:length(switchdt)
     switchdts(ii)=datenum(switchdt{ii},'yyyy-mm-dd HH:MM:SS');
end

for ii=1:length(bnds)
    bndsjs(ii)=datenum(bnds{ii},'yyyy-mm-dd HH:MM:SS');
end
%combine the fvs and plot.

for ii=1:length(pathvl)
    strcmd=['load ' pathvl{ii} catchvl{ii} '.mat']
    eval(strcmd);
    fvtmp=fvb;
    fvcomb=[fvcomb fvtmp];
end
 vals=getvals(fvcomb,2,'TRIG');
%now make a basic plot of all the fvs.

vlot=ceil((vals(:,1)-bndsjs(1)));
switchdays=reduce_dates(switchdts,bndsjs(1));
daysbnds=bndsjs(2), bndsjs(1))
figure

%plot(vlot, vals(:,2),'.')
%box off;
hold on;
%overlay mean and standard deviation
days=[1:daysbnds]

unique(vlot);

for ii=1:length(days)
    indtmp=find(vlot==days(ii))
    mnvl(ii)=mean(vals(indtmp,2))
    stdv(ii)=std(vals(indtmp,2))
%to add means and standard deviations.
    nnotes(ii)=length(indtmp)
    if (calc_prop_note)
        %set the initial file name, this is necessary to count total
        %segmented notes once....
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
plot(daystm,htrt,'o','Linewidth',5);
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








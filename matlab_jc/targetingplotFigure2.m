function [apctrace]=targetingplotFigure2(shortlong,Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff,notesizeshiftms)
    % shortlongfeeds in either Alldata(1).ind_shortnotes,
    % Alldata(1).ind_longnotes or Alldata(1).ind_allnotes
    
%%%%% OUTPUTS
    % slopes - linear slopes of decay over first 120pts (15ms on either side)
    % sdtoff - precision of targeting (1/std(targeting in ms))
    % mntoff - mean position of targeting
    % mnshift - position of maximal pitch shift
    


% prctile_cutoff=0.8;
% tmp_cutoff=24;
% shiftedmethod=1;
% notesizeshiftms=0;
notesizeshift=notesizeshiftms*8;

for i=1:length(Alldata)
    ons=Alldata(i).startnote;
    offs=Alldata(i).endnote;
    exp=Alldata(i).exp;
    shift_direction=Alldata(i).exp(1).direction;
    
    % get an indicator of the mean pitch in each 4-14hr group of song
    count=0;
    x=[];
    y=[];
    z=[];
    for j=1:length(exp)
        if exp(1).acsf(j)==1 && exp(1).baseline(j)==0
            count=count+1;
            x(count)=(exp(1).tfromwnonend(j)+exp(1).tfromwnonbegin(j))/2;
            y(count)=mean(mean(exp(j).selectedpitchcurves(ons+40:offs-40,:)'));
            z(count).toffs=exp(j).toffset;
            z(count).curves=exp(j).selectedpitchcurves;
        end
    end
    % 40 points is 5ms
    avgbaseline=mean(Alldata(i).baselineAC(ons+40:offs-40,:)');
    
    if strcmp(shift_direction,'up')
        up=1;
        for k=1:length(y)
            [extremum,ind]=max(y);
        end
    else
        up=0;
        for k=1:length(y)
            [extremum,ind]=min(y);
        end
    end
    maxdist=abs(extremum-avgbaseline);
    shiftedinds=[];
    if shiftedmethod==1;
    % prctile_cutoff method
        count1=0;
        for ll=1:length(y)
            if abs(y(ll)-avgbaseline)/maxdist > prctile_cutoff
                count1=count1+1;
                shiftedinds(count1)=ll;
            end
        end
    else
    % temporal proximity method
        count2=0;
        for kk=1:length(y)
            if abs(x(kk)-x(ind))<=tmp_cutoff
                count2=count2+1;
                shiftedinds(count2)=kk;
            end
        end
    end
    
    % get all toffs during shifting (before it's shifted)
    toffsets=[];
    for jj=1:length(y)
        if jj<=ind
            toffsets=[toffsets;z(jj).toffs];
        end
    end
    to(i).dat=toffsets;
    bigtarget(i)=(median(toffsets-16)-(ons+notesizeshift))/((offs-notesizeshift)-(ons+notesizeshift)); 

    % get all pitch curves during shifted 
    pcurves=[];
    for iii=1:length(shiftedinds)
        index=shiftedinds(iii);
        pcurves=[pcurves z(index).curves];
    end
    
    avgpcurvetrace(i).dat=(mean(pcurves(ons:offs,:)')-mean(Alldata(i).baselineAC(ons:offs,:)')); %./std(Alldata(i).baselineAC(ons:offs,:)');
    if up==1
        [a(i),bshift(i)]=max(avgpcurvetrace(i).dat);
    else
        [a(i),bshift(i)]=min(avgpcurvetrace(i).dat);
    end
    bigshift(i)=(bshift(i)-notesizeshift)/((offs-notesizeshift)-(ons+notesizeshift));
end

%figure;plot(bigtarget,bigshift,'*')

apctrace=zeros(length(Alldata),1000);
for i=1:length(Alldata)
    apctrace(i,1:length(avgpcurvetrace(i).dat))=(avgpcurvetrace(i).dat/a(i));
end

g=8;

abc=zeros(length(Alldata),1000);

for aaa=1:length(shortlong);
    n=shortlong(aaa);
left=bshift(n)-1;
right=length(avgpcurvetrace(n).dat)-bshift(n);
abc(n,1)=apctrace(n,bshift(n));
if right>left
    for m=2:left
        abc(n,m)=mean([apctrace(n,bshift(n)-m) apctrace(n,bshift(n)+m)]);
    end
    for m=left+1:right
        abc(n,m)=apctrace(n,bshift(n)+m);
    end
else
    for m=2:right
        abc(n,m)=mean([apctrace(n,bshift(n)-m) apctrace(n,bshift(n)+m)]);
    end
    for m=right+1:left
        abc(n,m)=apctrace(n,bshift(n)-m);
    end
end
end
    

for nn=1:size(abc,2)
    indices=find(abc(:,nn)~=0);
    avcurve(nn)=mean(abc(indices,nn));
    avup(nn)=mean(abc(indices,nn))+std(abc(indices,nn))/sqrt(length(shortlong));
    avdn(nn)=mean(abc(indices,nn))-std(abc(indices,nn))/sqrt(length(shortlong));
end
if length(shortlong)<10
    xaxis=(0:1/8:159/8);
    plot(xaxis,avcurve(1:160),'r');hold on;plot(xaxis,avup(1:160),'r');hold on;plot(xaxis,avdn(1:160),'r')
else
    xaxis=(0:1/8:319/8);
    plot(xaxis,avcurve(1:320),'k');hold on;plot(xaxis,avup(1:320),'k');hold on;plot(xaxis,avdn(1:320),'k')
end

xlim([0 85])

% Get all toffsets and graph on the same plot
toall=[];
for i=1:length(Alldata)
    toall=[toall;to(i).dat-mean(to(i).dat)];
end
histo=hist(abs(toall)/8,100);
%hold on;plot(histo/max(histo),'r')



% Do it for right over 1sigma
firstover1sig=[8 9 9 5 13 13 11 4 7 7 9 7 6 11 12 8];

to(1).dat=to(1).dat(1:120);
to(7).dat=to(7).dat(1:180);
to(9).dat=to(9).dat(400:600);
to(12).dat=to(12).dat(1:180);
to(14).dat=to(14).dat(1:180);
to(15).dat=to(15).dat(1:180);


x=(0:1/8:119/8);
for i=1:length(shortlong)
    n=shortlong(i);
    k=polyfit(x,abc(n,1:120),1);
    slopes(i)=-k(1);
    sdtoff(i)=1/std(to(n).dat);
    mntoff(i)=bigtarget(i);
    mnshift(i)=bigshift(i);
end
g=8;


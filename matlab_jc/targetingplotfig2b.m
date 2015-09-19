function [avgpcurvetrace]=targetingplotfig2b(Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff,notesizeshiftms)

% prctile_cutoff=0.8;
% tmp_cutoff=24;
% shiftedmethod=1;
% notesizeshiftms=0;
notesizeshift=notesizeshiftms*8;
firstover1sig=[8 9 9 5 13 13 11 4 7 7 9 7 6 11 12 8];


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
    
    % This takes into account the datasets that are longer than 6 hours.
    avgfirstpcurvetrace(i).dat=(mean(exp(firstover1sig(i)).selectedpitchcurves(ons:offs)')-mean(Alldata(i).baselineAC(ons:offs,:)'));
    if i==1    avgfirstpcurvetrace(i).dat=(mean(exp(firstover1sig(i)).selectedpitchcurves(ons:offs,1:20)')-mean(Alldata(i).baselineAC(ons:offs,:)'));
    else if i==6     avgfirstpcurvetrace(i).dat=(mean(exp(firstover1sig(i)).selectedpitchcurves(ons:offs,100:183)')-mean(Alldata(i).baselineAC(ons:offs,:)'));
        else if i==12        avgfirstpcurvetrace(i).dat=(mean(exp(firstover1sig(i)).selectedpitchcurves(ons:offs,1:25)')-mean(Alldata(i).baselineAC(ons:offs,:)'));
            else if i==16       avgfirstpcurvetrace(i).dat=(mean(exp(firstover1sig(i)).selectedpitchcurves(ons:offs,1:40)')-mean(Alldata(i).baselineAC(ons:offs,:)'));
                end
            end
        end
    end
    
    
    if up==1
        [a(i),bshift(i)]=max(avgfirstpcurvetrace(i).dat);
    else
        [a(i),bshift(i)]=min(avgfirstpcurvetrace(i).dat);
    end
    bigshift(i)=(bshift(i)-notesizeshift)/((offs-notesizeshift)-(ons+notesizeshift));
end

%figure;plot(bigtarget,bigshift,'*')

apctrace=zeros(length(Alldata),1000);
for i=1:length(Alldata)
    apctrace(i,1:length(avgfirstpcurvetrace(i).dat))=(avgfirstpcurvetrace(i).dat/a(i));
end

g=8;

abc=zeros(length(Alldata),1000);

for aaa=1:length(Alldata(1).ind_longnotes);   % JUST THE LONG NOTES
    n=Alldata(1).ind_longnotes(aaa);
left=bshift(n)-1;
right=length(avgfirstpcurvetrace(n).dat)-bshift(n);
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
    avup(nn)=mean(abc(indices,nn))+std(abc(indices,nn));
    avdn(nn)=mean(abc(indices,nn))-std(abc(indices,nn));
    xaxis(nn)=nn/8;
end
figure;plot(xaxis,avcurve);hold on;plot(xaxis,avup);hold on;plot(xaxis,avdn)
xlim([0 85])

% Get all toffsets and graph on the same plot
toall=[];
for ii=1:length(Alldata(1).ind_longnotes)
    i=Alldata(1).ind_longnotes(ii);
    toall=[toall;to(i).dat(1:100)-mean(to(i).dat(1:100))];
end
histo=hist(abs(toall)/8,100);
hold on;plot(histo/max(histo),'r')

% Do it for right over 1sigma



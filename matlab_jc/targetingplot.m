function [bigtarget,bigshift]=targetingplot(Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff,notesizeshiftms)

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
        [a,bigshift(i)]=max(avgpcurvetrace(i).dat);
    else
        [a,bigshift(i)]=min(avgpcurvetrace(i).dat);
    end
    bigshift(i)=(bigshift(i)-notesizeshift)/((offs-notesizeshift)-(ons+notesizeshift));
end

plot(bigtarget,bigshift,'*','Color','r')



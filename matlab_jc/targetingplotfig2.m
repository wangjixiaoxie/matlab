function [avgpcurvetrace]=targetingplotfig2(Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff,notesizeshiftms)
% Uses both long and short stack notes and returns both avgpcurvetrace and
% shif/sd which show the weak relationship between std of targeting and sharpness of
% pitch peak.

% prctile_cutoff=0.8;
% tmp_cutoff=24;
% shiftedmethod=1;
% notesizeshiftms=0;
notesizeshift=notesizeshiftms*8;

for i=1:length(Alldata)
    % Information for plotting things as proportion of the note's length.
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

for i=1:length(Alldata)
    abc(i,:)=onesided(Alldata(i).);
end
g=7;
% for nn=1:size(abc,2)
%     indices=find(abc(:,nn)~=0);
%     avcurve(nn)=mean(abc(indices,nn));
%     avup(nn)=mean(abc(indices,nn))+std(abc(indices,nn));
%     avdn(nn)=mean(abc(indices,nn))-std(abc(indices,nn));
%     xaxis(nn)=nn/8;
% end
% plot(xaxis,avcurve);hold on;plot(xaxis,avup);hold on;plot(xaxis,avdn)
% xlim([0 85])
% 
% % Get all toffsets and graph on the same plot
% toall=[];
% for i=1:length(Alldata)
%     toall=[toall;to(i).dat-mean(to(i).dat)];
% end
% histo=hist(abs(toall)/8,100);
% %hold on;plot(histo/max(histo),'r')
% 
% 
% 
% % % Calculate 
% for j=1:160 % 5ms to 20ms
%     for i=1:length(Alldata)
%         shif(i)=1-abc(i,j); % amount of peak reduced by 10ms on either side
%         sd(i)=std(to(i).dat);
%     end
%     gg=corrcoef(sd,shif);
%     ccoef(j)=gg(2);
% end

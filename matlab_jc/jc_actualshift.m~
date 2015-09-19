function [actshift,normshift,pcnormshift,maxshiftlocus]=jc_actualshift(Alldata,shiftedmethod,prctile_cutoff,tmp_cutoff)
% Returns the actual shift, normalized shift, and peak-centered normalized
% shift for a single experiment.  maxshiftlocus is the location of maximal
% shift as a proportion of the length of the note.

% jc_actualshift(Alldata(6),1,0.8,24);

% Information for plotting things as proportion of the note's length.
ons=Alldata.startnote;
offs=Alldata.endnote;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp=Alldata.exp;
shift_direction=Alldata.exp(1).direction;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get an metric of the mean pitch in each 4-14hr group of song
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
% Take the mean of the middle of the baseline notes --- 40 points is 5ms
avgbaseline=mean(Alldata.baselineAC(ons+40:offs-40,:)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%% Determines when note is "shifted" %%%%%%%%%%%%%
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

 
%%%%%%%%%% Gets mean of pitch curves during shifted%%%%%%%%
    pcurves=[];
    for iii=1:length(shiftedinds)
        index=shiftedinds(iii);
        pcurves=[pcurves z(index).curves];
    end

    actshift=(mean(pcurves(ons:offs,:)')-mean(Alldata.baselineAC(ons:offs,:)')); 

%%%%%%%%%%% Calculates location of maximum shift as proportion of note
    if up==1
        [a,bshift]=max(actshift);
    else
        [a,bshift]=min(actshift);
    end
    maxshiftlocus=(bshift)/((offs)-(ons));

%%%%%%%%%%% Normalizes by direction and magnitude
    normshift=zeros(1,1000);
    normshift(1:length(actshift))=(actshift/a);
%%%%%%%%%%% Normalizes by side
    pcnormshift=onesided(normshift,ons,offs);

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

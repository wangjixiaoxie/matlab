%4.16.09, purpose is to calculate mean, stderror of matrix, excluding zero
function [outmean, outsderror,pvl]=calcmeanstder2(matin)

if(isempty(matin))
    outmean=[];
    outsderror=[];
else
%initialize to zero
for ii=1:length(matin(1,:))
    outmean(ii)=0;
    outsderror(ii)=0;
end

for ii=1:length(matin(1,:))
    notnan=find((matin(:,ii)~=0)&(~isnan(matin(:,ii))));
    if(~isempty(notnan))
        outmean(ii)=mean(matin(notnan,ii));
%         [h,p]=ttest(matin(notnan,ii));
%         pvl(ii)=p;
pvl=1;
    end
end
%stderror I need to loop through column by column

for ii=1:length(matin(1,:))
    notnan=find((matin(:,ii)~=0)&(~isnan(matin(:,ii))));
    stdvcur=std(matin(notnan,ii));
    stdercur=stdvcur/sqrt(length(notnan));
    if(~isempty(notnan))
    outsderror(ii)=stdercur;
    end
end
end
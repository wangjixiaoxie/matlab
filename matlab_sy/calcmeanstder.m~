%4.16.09, purpose is to calculate mean, stderror of matrix, excluding nan
function [outmean, outsderror]=calcmeanstder(matin)

notnan=~isnan(matin);
matin(~notnan)=0;
howmany=sum(notnan);

%for mean
columnTot=sum(matin);
outmean=columnTot./howmany;

%stderror I need to loop through column by column

for ii=1:length(matin(1,:))
    notnan=~isnan(matin(:,ii));
    stdvcur=std(matin(notnan,ii));
    stdercur=stdvcur/sqrt(length(notnan));
    outsderror(ii)=stdercur;
end
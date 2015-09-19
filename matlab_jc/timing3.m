function [timevals]=timing3(fvals)
% fvals=findwnoteJC('batchnotes','a','','',0,[3000 4200],8000,1,'obs0',1);
% mk_tempf('batchfiles',templaJC1,2,'obs0');
% get_trigt2('batchfiles',cntrngJC1,0.35,128,1,1);
% vals=triglabel('batchfiles','a',1,1,0,1);
% **** for leap years, change monthdays
% 
if isempty(fvals)
    timevals=[];
    return;
end
monthdays=[31 28 31 30 31 30 31 31 30 31 30 31];
monthsums(1)=0;
for i=2:length(monthdays)
    monthsums(i)=sum(monthdays(1:i-1));
end
count=0;
fvalscntr=0;

% For each song
for i=1:length(fvals)
    hourplc=findstr(fvals(i).fn,'_');
    otherplc=findstr(fvals(i).fn,'.');
    if otherplc(1)-hourplc(2)<6 % evtaf old
        
            timevals(i)=24*(monthsums(fvals(i).month)+(fvals(i).day-1))+fvals(i).hour;
    else  % Evtaf4
            thismonth=str2num(fvals(i).fn(hourplc(1)+3:hourplc(1)+4));
            thisday=str2num(fvals(i).fn(hourplc(1)+1:hourplc(1)+2));
            hours=str2num(fvals(i).fn(hourplc(2)+1:hourplc(2)+2))+fvals(i).hour/60;
            timevals(i)=24*(monthsums(thismonth)+(thisday-1))+hours;
    end
end

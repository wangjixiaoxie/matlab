function [hour,day,month,year]=fn2datesec(ebinfn);
%[hour,day,month,year]=fn2date(ebinfn);
%

p = findstr(ebinfn,'.');
p = p(end-1);
hr = ebinfn([(p(end)-6):(p(end)-3)]);
dt = ebinfn([(p(end)-13):(p(end)-8)]);

hour = str2num(hr(1:2)) + str2num(hr(3:4))/60.0;
day   = str2num(dt(1:2));
month = str2num(dt(3:4));
year  = 2000 + str2num(dt(5:6));
return;

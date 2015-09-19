function [hour,day,month,year]=fn2date(ebinfn);
%[hour,day,month,year]=fn2date(ebinfn);
%

p = findstr(ebinfn,'.');
p = p(end-1);
hr = ebinfn([(p(end)-4):(p(end)-1)]);
dt = ebinfn([(p(end)-11):(p(end)-6)]);

hour = str2num(hr(1:2)) + str2num(hr(3:4))/60.0;
day   = str2num(dt(1:2));
month = str2num(dt(3:4));
year  = 2000 + str2num(dt(5:6));
return;

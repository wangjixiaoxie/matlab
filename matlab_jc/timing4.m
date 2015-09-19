function [timevals]=timing4(fvals)
% fvals=findwnoteJC('batchnotes','a','','',0,[3000 4200],8000,1,'obs0',1);
% mk_tempf('batchfiles',templaJC1,2,'obs0');
% get_trigt2('batchfiles',cntrngJC1,0.35,128,1,1);
% vals=triglabel('batchfiles','a',1,1,0,1);
% **** for leap years, change monthdays
% 
mv=[0 31 59 90 120 151 181 212 243 273 304 334];
% For each song
for i=1:length(fvals)
    a=find(fvals(i).fn=='_');
    monthval=str2num(fvals(i).fn(a(1)+3:a(1)+4));
    months=mv(monthval)*24;
    days=str2num(fvals(i).fn(a(1)+1:a(1)+2))*24;
    hours=str2num(fvals(i).fn(a(2)+1:a(2)+2));
    minutes=str2num(fvals(i).fn(a(2)+3:a(2)+4))/60;
    if strcmpi(fvals(i).fn(a(2)+5:a(2)+6),'.-')
        seconds=str2num([fvals(i).fn(a(2)+5) fvals(i).fn(a(2)+7)])/3600;
    else
        seconds=str2num(fvals(i).fn(a(2)+5:a(2)+6))/3600;
    end
    timevals(i)=months+days+hours+minutes+seconds;
end

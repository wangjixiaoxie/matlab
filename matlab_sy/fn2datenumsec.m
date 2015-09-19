function [dtnm]=fn2date(fn);
%dtnm=fn2date(fn);
%

ISWAV=0;
p = findstr(fn,'.cbin');
if (length(p)<1)
    p=findstr(fn,'.ebin');
end
if (length(p)<1)
    p=findstr(fn,'.rec');
end
if (length(p)<1)
    p=findstr(fn,'.bbin');
end
if (length(p)<1)
    p=findstr(fn,'.wav');
    ISWAV=1;
end

if (length(p)<1)
    disp(['not rec cbin, ebin or bbin file?']);
    return;
end

if (ISWAV==0)
    p2=findstr(fn,'.');
    p3=find(p2==p(end));
    if (length(p3)<1)|(p3(1)<2)|(length(p2)<2)
        disp(['weird fn = ',fn]);
        return;
    end
    p = p2(p3-1);
    hr   = fn([(p(end)-6):(p(end)-5)]);
    mnut = fn([(p(end)-4):(p(end)-3)]);
    dy   = fn([p(end)-13:p(end)-12]);
    mnth = fn([p(end)-11:p(end)-10]);
    yr   = fn([p(end)-9:p(end)-8]);
    dt   = fn([(p(end)-13):(p(end)-8)]);
    scnd='0';
else
    pp=findstr(fn,'_');
    mnth=fn(pp(2)+1:pp(3)-1);
    dy =fn(pp(3)+1:pp(4)-1);
    hr=fn(pp(4)+1:pp(5)-1);
    mnut=fn(pp(5)+1:pp(6)-1);
    ppp=findstr(fn,'.wav');
    scnd=fn(pp(6)+1:ppp(1)-1);
    ff=dir(fn);
    yr=ff(1).date;
    ppp=findstr(yr,'-');
    yr = yr(ppp(2)+1:ppp(2)+4);
    
    dt=[dy,mnth,yr(3:4)];
end

hour = str2num(hr) + str2num(mnut)/60.0 +str2num(scnd)/3600.0;
day   = str2num(dy);
month = str2num(mnth);
year  = 2000 + str2num(yr);
dtnm = datenum([year,month,day,str2num(hr),str2num(mnut),str2num(scnd)]);
return;

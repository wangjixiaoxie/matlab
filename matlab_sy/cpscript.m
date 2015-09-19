function [sm,sp,t,f]=cpscript(destdir,birdnum,datenum)

timebnds={'07' '08' '09' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21'}

for ii=1:length(datenum)
    curdt=datenum(ii);
    for jj=1:length(timebnds)
        curtm=timebnds{jj};
        if (curdt>9)
            cmd=['!cp -u *' num2str(birdnum) '_' num2str(curdt) '*_' num2str(curtm) '* ' destdir]
        else
             cmd=['!cp -u *' num2str(birdnum) '_0' num2str(curdt) '*_' num2str(curtm) '* ' destdir]
        end
             eval(cmd);   
    end
end
function tvals=timing(vals,fvals)
% fvals=findwnoteJC('batchnotes','a','','',0,[3000 4200],8000,1,'obs0',1);
% mk_tempf('batchfiles',templaJC1,2,'obs0');
% get_trigt2('batchfiles',cntrngJC1,0.35,128,1,1);
% vals=triglabel('batchfiles','a',1,1,0,1);
% **** for leap years, change monthdays
% 
monthdays=[31 28 31 30 31 30 31 31 30 31 30 31];
monthsums(1)=0;
for i=2:length(monthdays)
    monthsums(i)=sum(monthdays(1:i-1));
end
count=0;
fvalscntr=0;

% For each song
for i=1:size(vals,1)
    numfvs=vals(i,2);
    fvalscntr=fvalscntr+numfvs;
    if fvalscntr>0
        songtime(i)=24*(monthsums(fvals(fvalscntr).month)+(fvals(fvalscntr).day-1))+fvals(fvalscntr).hour;
    else 
        vals(i,1)=0;
    end
    if vals(i,1)>0
        for j=1:vals(i,1)
            tvals(count+j)=songtime(i);
        end
        count=count+j;
    end

end
g=7;
    
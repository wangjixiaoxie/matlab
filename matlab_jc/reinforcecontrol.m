function [postneg,postescape]=reinforcecontrol(meanpitch)
a=1;
b=1;
c=1;
mp=meanpitch;
for i=3:length(mp)
    if(mp(i-1)>2380)
        postneg(a)=mp(i);
        a=a+1;
    else 
        postescape(b)=mp(i);
        b=b+1;
    end
end

    
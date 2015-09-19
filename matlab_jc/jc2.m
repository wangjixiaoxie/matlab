function [gp,np]=jc2(p714a)
for i=1:length(p714a)
    if std(p714a(i).pitches(750:850))<80
        gp(i)=median(p714a(i).pitches(320:400));
        np(i)=0;
    else
        gp(i)=0;
        np(i)=median(p714a(i).pitches(320:400));
    end
end
function [postneg,postescape,postex]=reinforce(boomed,meanpitch)
a=1;
b=1;
c=1;
mp=meanpitch;
for i=3:length(mp)
    if boomed(i)<100
    if(boomed(i-1)>100)
        postneg(a)=mp(i);
        a=a+1;
    else if(boomed(i-2)>100)
            postescape(b)=mp(i);
            b=b+1;
        else postex(c)=mp(i);
            c=c+1;
        end
    end
    end
end

    
function [diff,zeroe]=zeros0(xcresids)
k=1;

for i=1:length(xcresids)-1
    if xcresids(i)==0
        zeroe(k)=i;
        k=k+1;
    else if xcresids(i)<0 && xcresids(i+1)>0 || xcresids(i)>0 && xcresids(i+1)<0
        zeroe(k)=i+0.5;
        k=k+1;
        end
    end
end
for i=1:length(zeroe)
    if zeroe(i)<700
        a(1)=zeroe(i);
    end
    if zeroe(i)>700
        a(2)=zeroe(i);
        break
    end
end
diff=a(2)-a(1);



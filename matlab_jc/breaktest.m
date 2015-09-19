function b=breaktest
a=[1 2 3 4 5 6];
for i=1:length(a)
    if a(i)<3.5
        b(1)=a(i);
    end
    if a(i)>3.5
        b(2)=a(i);
        break
    end
end
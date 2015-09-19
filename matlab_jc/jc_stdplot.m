function thresh=jc_stdplot(input,width)
long=length(input)-width;
for i=1:long
    output(i)=std(input(i:i+width));
    if output(i)<50
        thresh=i;
        break
    end
end
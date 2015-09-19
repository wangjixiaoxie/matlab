function thresh=jc_getthresh(input,width)

%This can be used to align the notes by using a decrease in the 
%standard deviation below some threshold as an indicator of note onset.
long=length(input)-width;
for i=1:long
    output(i)=std(input(i:i+width));
    if output(i)<25
        thresh=i;
        break
    end
end
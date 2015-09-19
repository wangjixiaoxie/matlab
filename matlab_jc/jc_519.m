function avgnote=jc_519(array)
%Input is the output array from findwnote2

%This loop finds the average note.
for i=1:length(array(1).datt)
    sumnote(i)=0;
    for j=1:length(array)
        sumnote(i)=sumnote(i)+array(j).datt(i);
    end
    avgnote(i)=sumnote(i)/length(array);
end


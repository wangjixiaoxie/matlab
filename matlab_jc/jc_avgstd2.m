function jc_avgstd(arrayfile,width)

%Get average note
for i=1:length(arrayfile(1).pitches)
    sumnote(i)=0;
    for j=1:length(arrayfile)
        sumnote(i)=sumnote(i)+arrayfile(j).pitches(i);
    end
    avgnote(i)=sumnote(i)/length(arrayfile);
end

%This can be used to align the notes by using a decrease in the 
%standard deviation below some threshold as an indicator of note onset.
long=length(arrayfile(1).pitches)-width-1;

for i=1:long
    stand(i)=0;
    for j=1:length(arrayfile)
        stand(i)=stand(i)+std(arrayfile(j).pitches(i:i+width));
    end
    stand(i)=stand(i)/length(arrayfile);
end
figure; plot(stand)
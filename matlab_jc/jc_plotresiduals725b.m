function normalized=jc_plotresiduals725b(pitchdata)



%Get the residuals
for i=1:length(pitchdata(1).pitches)
    t(i)=0;
    for j=1:length(pitchdata)
        matr(i,j)=pitchdata(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
end

for i=1:length(pitchdata)
    for j=1:length(pitchdata(1).pitches)
        normalized(j,i)=matr(j,i)/averaged(j);
        normalized(j,i)=normalized(j,i)-1;
    end
end



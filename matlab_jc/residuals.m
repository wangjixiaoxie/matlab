function normalized=residuals(pitchmatrix)
pitchdata=pitchmatrix;
%Get the residuals
for i=1:size(pitchdata,1);avg(i)=mean(pitchdata(i,:));end

for i=1:size(pitchdata,2)
    for j=1:size(pitchdata,1)
        normalized(j,i)=pitchdata(j,i)/avg(j);
        normalized(j,i)=normalized(j,i)-1;
    end
end



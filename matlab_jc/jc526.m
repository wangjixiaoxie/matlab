function deviants=jc526(input)
for i=1:30
    meanest=mean(input(i,:));
    for j=1:3
        deviants(i,j)=input(i,j)-meanest;
    end
end
function [bestmin,marker,combo3]=PCAclassify3(SCORE,size1)
% for the first 10 components, determine the goodness of classification
% What variation is present in UD but not D?
numsteps=400;

for i=1:10
    mmax=max(abs(SCORE(:,i)));
    stepwidth=(mmax)/numsteps;
    for ii=1:numsteps
        divider=0+stepwidth*ii;
        errors(ii)=length(find(abs(SCORE(1:size1,i))<divider))+length(find(abs(SCORE(size1+1:end,i))>divider));
    end
    [bestmin(i),mark]=min(errors);
    marker(i)=mark*stepwidth;
end

for j=2:10
    for k=2:10
        for kk=2:10
        combo3(j-1,k-1,kk-1)=length(find(abs(SCORE(find(abs(SCORE(find(abs(SCORE(1:size1,j))>marker(j)),k))>marker(k)),kk)>marker(kk))))+length(find(abs(SCORE(find(abs(SCORE(find(abs(SCORE(size1+1:end,j))<marker(j)),k))<marker(k)),kk)<marker(kk))));
        end
    end
end
1-(min(min(min(combo3)))/size(SCORE,1))

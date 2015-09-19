function mini=PCAclassify2(SCORE,size1)
% for the first 10 components, determine the goodness of classification
% What variation is present in UD but not D?
for i=1:10
    mmin=min(SCORE(:,i));
    mmax=max(SCORE(:,i));
    stepwidth=(mmax-mmin)/100;
    % choose random boundaries
    for ii=1:1000
        rb1(ii)=rand*100*stepwidth+mmin;
        rb2(ii)=rand*100*stepwidth+mmin;
        randbound1=rb1(ii);
        randbound2=rb2(ii);
        if randbound1>randbound2
            errors(ii)=length([find(SCORE(find(SCORE(1:size1,i)>randbound2),i)<randbound1);(find(SCORE(size1+1:end,i)<randbound2));(find(SCORE(size1+1:end,i)>randbound1))]);
        else
            errors(ii)=length(find(SCORE(find(SCORE(1:size1,i)<randbound2),i)>randbound1))+length(find(SCORE(size1+1:end,i)>randbound2))+length((find(SCORE(size1+1:end,i)<randbound1)));
        end
    end
    [mini(i),marker]=min(errors);
    
end
        
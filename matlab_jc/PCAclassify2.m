function [bestmin,marker,combo2]=PCAclassify2(SCORE,size1)
% for the first 10 components, determine the goodness of classification
% What variation is present in UD but not D?
numsteps=400;

% for i=1:10
%     mmax=max(abs(SCORE(:,i)));
%     stepwidth=(mmax)/numsteps;
%     for ii=1:numsteps
%         divider=0+stepwidth*ii;
%         errors(ii)=length(find(abs(SCORE(1:size1,i))<divider))+length(find(abs(SCORE(size1+1:end,i))>divider));
%     end
%     [bestmin(i),mark]=min(errors);
%     marker(i)=mark*stepwidth;
% end

counter=0;
for jj=2:10
for j=2:10
    for k=2:10
        for kk=2:10
            counter=counter+1;
            sumX=abs(SCORE(:,jj))+abs(SCORE(:,j))+abs(SCORE(:,k))+abs(SCORE(:,kk));    
            mmax=max(sumX);
            stepwidth=(mmax)/numsteps;
        for ii=1:numsteps
            divider=0+stepwidth*ii;
            errors(ii)=length(find(abs(sumX(1:size1))<divider))+length(find(abs(sumX(size1+1:end))>divider));
        end
        g(counter)=min(errors);
        end
    end
end
end
%%%%%%%%%%%
%%%%%%%%%%%
%
numsteps=400;
mmax=5;
stepwidth=mmax/numsteps;
for ii=1:numsteps
    divider=0+stepwidth*ii;
    errors(ii)=length(find(mm(1:UDend)<divider))+length(find(mm(UDend+1:end)>divider));
end
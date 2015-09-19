function ravg=db_runningaverage(input,window,dim)
%input a matrix, window size, and the dimension, and it will give you a
%running average

if (~exist('dim'))
    dim = 1;
end
    
if dim == 1
    for i=1:length(input)-window
        ravg(i,:)=mean(input(i:i+window,:),dim);
    end
elseif dim ==2
    for i=1:length(input)-window
        ravg(:,i)=mean(input(:,i:i+window),dim);
    end
end


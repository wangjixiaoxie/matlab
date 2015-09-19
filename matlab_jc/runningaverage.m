function ravg=runningaverage(input,window,dim)

if (~exist('dim'))
    dim = 1;
end
    

for i=1:length(input)-window
    ravg(i)=mean(input(i:i+window),dim);
end


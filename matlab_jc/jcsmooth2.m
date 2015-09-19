function ravg=runningaverage(input,window)
for i=1:length(input)-window
    ravg(i)=mean(input(i:i+window));
end
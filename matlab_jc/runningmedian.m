function ravg=runningmedian(input,window)
for i=1:length(input)-window
    ravg(i)=median(input(i:i+window));
end
function ravg=runningprctile(input,prct,window)
for i=1:length(input)-window
    ravg(i)=prctile(input(i:i+window),prct);
end
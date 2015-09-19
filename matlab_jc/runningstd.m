function ravg=runningstd(input,window)
for i=1:length(input)-window
    ravg(i)=std(input(i:i+window));
end
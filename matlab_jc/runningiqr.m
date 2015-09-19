function IQR=runningiqr(input,window)
for i=1:length(input)-window
    IQ1=prctile(input(i:i+window),25);
    IQ3=prctile(input(i:i+window),75);
    IQR(i)=(IQ3-IQ1);
end
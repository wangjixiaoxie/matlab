function y=jcfiltfiltbatch(norm)
for i=1:size(norm,2)
    x=norm(:,i);
    [b,a] = butter(4,300/1000,'high');
    y(:,i)=filtfilt(b,a,x);
end
    
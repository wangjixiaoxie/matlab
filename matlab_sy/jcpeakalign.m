function middles=jcpeakalign(fvals,first,last,width)

fs=32000;
%Smooth the raw data
for i=1:length(fvals)
    [holder]=SmoothData(fvals(i).datt,fs,1);
    smooth(i,:)=log(holder);
end

for i=1:25
    for j=first+width:last
        ss(i,j)=sum(smooth(i,j-width:j));
    end
    [a,b]=max(ss(i,:));
    middles(i)=b-width/2;
end
figure;hold on;
for i=1:25
    plot(jc_hilbert(syllable(i,middles(i)-999:middles(i)+9500),2000,5000,2))
end

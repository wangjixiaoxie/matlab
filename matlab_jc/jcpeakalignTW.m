function curves=jcpeakalignTW(fvals,first,last,width)

fs=32000;
%Smooth the raw data
for i=1:length(fvals)
    [holder]=SmoothData(fvals(i).datt,fs,1);
    smooth(i,:)=log(holder);
end

for i=1:50
    for j=first+width:last
        ss(i,j)=sum(smooth(i,j-width:j));
    end
    [a,b]=max(ss(i,:));
    middles(i)=b-width/2;
end
%figure;hold on;
for i=1:50
    curves(:,i)=(jc_hilbert(smooth(i,middles(i)-699:middles(i)+2500),2600,4200,1));
end

function middles=jcpeakalign(fvals,first,last,width)

fs=32000;
%Smooth the raw data
for i=1:length(fvals)
    [holder]=SmoothData(fvals(i).datt,fs,1);
    smooth(i,:)=log(holder);
end
first=[500 4000 5000];
last=[3000 5000 6000];
middle=1000;
for k=1:3
    clear ss
    for i=1:25
        for j=first(k)+width:last(k)
            ss(i,j)=sum(smooth(i,j-width:j));
        end
        [a,b]=max(ss(i,:));
        middles(i)=b-width/2;
    end
    for i=1:25
        note(k).data(i,:)=jc_hilbert(smooth(i,middles(i)-899:middles(i)+2500),2000,5000,2);
    end
end
notes=[note(1).data(:,850:1850),note(2).data(:,850:1850),note(3).data(:,1000:2000)];
figure;plot(notes')
%%%% DTWarping
x1=B(1).selectedpitchcurves(450:749,2)-mean(B(1).selectedpitchcurves(450:749,1:19)')';
x2=B(2).selectedpitchcurves(450:749,1)-mean(B(2).selectedpitchcurves(450:749,1:19)')';
y1=fft(x1);
y2=fft(x2);
m1=abs(y1);
m2=abs(y2);
p1=unwrap(angle(y1));
p2=unwrap(angle(y2));
f=(0:length(y1)-1)'*100/length(y1);
figure;plot(f,m1)
hold on;plot(f,m2,'r')
figure;plot(f,p1*180/pi)
hold on;plot(f,p2*180/pi,'r')

pts=zeros(29,10);
for j=1:29
asw=answer(:,j);
for i=2:length(asw)
    dt(i)=asw(i)-asw(i-1);
end

count=0;
for i=1000:length(dt)
    if ((dt(i)>0 && dt(i-1)<0) || (dt(i)<0 && dt(i-1)>0))
        count=count+1;
        pts(j,count)=i;
    end
end
end



answer16=jc_hilbert(shifted2(:,1:15000),2000,4000,16);

figure;hold on;
gg=mean(answer16');
for i=1:size(answer16,2)
    chunk=answer16(pts(i,1):pts(i,2),i);
    long=length(chunk);
    y=resample(chunk,2000,long);
    points1(i,:)=(y-gg(pts(i,1):pts(i,1)+1999)');
end

for i=1:29
plot(answer(pts(i,1):pts(i,2),i))
end
for i=1:29
[a,b]=max(answer16(1000:2000,i));
c(i)=b+1000;
end
clear chunk
for i=1:size(answer16,2)
    chunk(:,i)=answer16(c(i):c(i)+10000,i);
end
for i=1:size(answer16,2)
    curves(:,i)=chunk(:,i)'-mean(chunk');
end
%%%%%% STEP 1: filter the original data
lowp=2800; % below band we care about power in - this defines note borders
highp=3800; % above band we care about power in - this defines note borders
samplerate=32000;
numnotes=size(RepeatBirds.shifted(2).Allrawdata,1);
% filter raw data and do hilbert transform
[b,a]=butter(4,[lowp/(0.5*samplerate) highp/(0.5*samplerate)],'bandpass');
for i=1:numnotes
    filter=filtfilt(b,a,RepeatBirds.shifted(2).Allrawdata(i,:));
    hilb4(i,:)=jc_hilbert(filter,lowp,highp,4);
    filtered(i,:)=filter;
end
%%%%% STEP 2: determine note onset and offset based on 50% of max raw trace
[smoothed]=RepeatSmooth


%%%%%%%%   

syllable=RepeatBirds.shifted(2).Allrawdata(1:10,:);
sonogram2=sono(syllable,1024,1020,2,2100,2800,[1 2 3],'obs0',1);
j=1;
for i=1:5000
    ste(i,j)=sum(log2(sonogram2(j).data(:,i)));
end

%%%%%
for i=1:10
    first=4000;
    last=7500;
    for j=first+1000:last
        ss(i,j)=sum(smoothed(i,j-1000:j));
    end
    [a,b(i)]=max(ss(i,:));
end
figure;hold on;
for i=1:10
    plot(jc_hilbert(syllable(i,b(i)-1599:b(i)+500),2000,5000,2))
end
    


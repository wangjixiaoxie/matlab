% SequenceTrack
tvalsWNON914dom=timing4(fvalsWNON914dom);
tvdom=tvalsWNON914dom;
tvrec=tvalsWNON914rec;
sizedom=size(tvdom,2);
tvall=[tvalsWNON914dom tvalsWNON914rec];
[b,ind]=sort(tvall);
width=50;
clear probdom0
for i=1:floor(length(tvall)/width)
    numdom=length(find(ind((i-1)*(width)+1:i*width)<sizedom));
    probdom0(i)=numdom/width;
end
figure;plot(probdom0)
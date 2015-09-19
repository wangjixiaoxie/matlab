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

%%%%% 1215 %%%%%%%%
x1=ACSF(3).residHILB2(1900:2500,21:50);
x2=INA(3).residHILB2(2000:2600,1:30);
for i=1:30
    a=ACSF(3).residHILB2(1900:2500,20+i);
    b=INA(3).residHILB2(2000:2600,0+i);
    z1=fft(a);
    z2=fft(b);
    mag1(i,:)=abs(z1);
    mag2(i,:)=abs(z2);
end
mag=mag1-mag2;
m=mean(mag)';
p=angle(z1);
re=m.*cos(p);     
im=m.*sin(p);
dat=complex(re,im);
figure;plot(real(ifft(dat)))

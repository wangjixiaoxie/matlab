function [xax,kk]=jcsd(pitch)
step=10;
width=50;
kk=zeros(1,200);
for j=1:10
for i=1:200
    start=(i-1)*step+1;
    kk(i)=kk(i)+std(pitch(j).pitches(start:start+width));
end
end
for i=1:200
    xax(i)=width/2+(i-1)*step;
end
function [outputdata]= jc_chirpy
for i=1:10000
    y(i)=sin(2*pi*i*(1+0.00001*i));
end
outputdata=y;
figure; plot(y)
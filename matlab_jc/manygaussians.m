function gauss_stor=manygaussians
a=0;
for i=1:10
    center=500-(i-10)*10;
    for j=0:4
        std=j*100+50;
        for k=1:4
            a=a+1;
            amp=0.001*(2*k+10);
            gauss_stor(a).center=center;
            gauss_stor(a).std=std;
            gauss_stor(a).amp=amp;
            gauss_stor(a).gauss=gaussian(amp,center,std);
        end
    end
end

            
function sumgauss=manysums
gauss=manygaussians;
a=0;
    for j=1:length(gauss)
        first=gauss(j).gauss;
        for k=1:length(gauss)
            a=a+1;
            second=gauss(k).gauss;
            sumgauss(a).gauss=first+second;
        end
    end
end

    
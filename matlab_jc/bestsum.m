function x=bestsum(y)
sumgauss=manysums;
diff=10^100;
place=911;
for i=1:length(sumgauss)
    difference=0;
    for j=1:length(y)
        difference=difference+(sumgauss(i).gauss(350+j)-y(j))^2;
    end
    if difference<diff
        place=i;
        diff=difference;
    end
end
x=sumgauss(place).gauss;
function F=jc_4gauss(x,xdata)
numgauss=4;
F=zeros(1,length(xdata));
    for j=1:numgauss
        first=(1+(j-1)*numgauss);
        second=(2+(j-1)*numgauss);
        third=(3+(j-1)*numgauss);
        for i=1:length(F)
            F(i)=F(i)+x(first)*exp(-(xdata(i)-x(second))^2/(2*x(third)^2));
        end
    end

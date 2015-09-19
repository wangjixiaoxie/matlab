function F=jc_SETgauss(x,xdata)
numgauss=3;
sig=150;
F=zeros(1,length(xdata));
    for j=1:numgauss
        first=(1+(j-1)*2);
        second=(2+(j-1)*2);
        %third=(3+(j-1)*3);
        for i=1:length(F)
            F(i)=F(i)+x(first)*exp(-(xdata(i)-x(second))^2/(2*sig^2));
        end
    end

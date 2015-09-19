function y=expfitjc(coef,x)
a=coef(1);
b=coef(2);
c=coef(3);

y=a*exp(-b*x)+c;
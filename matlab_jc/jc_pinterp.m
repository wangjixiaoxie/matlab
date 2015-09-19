function [output]=jc_pinterp(input)
a = 1;
b = [1/5 1/5 1/5 1/5 1/5];
output=filter(b,a,input);
   
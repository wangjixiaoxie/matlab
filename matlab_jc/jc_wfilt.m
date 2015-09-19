function [output]=jc_wfilt(input)

%This uses a basic moving average filter over 5 points
a = 1;
b = [1/5 1/5 1/5 1/5 1/5];
output=filter(b,a,input);
   
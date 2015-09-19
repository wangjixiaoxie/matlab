function [autocorr]=jc_xcorrnorm(residuals)
autocorr=[];
for i=1:size(residuals,2)
    autocorr(:,i)=xcorr(residuals(:,i))./max(xcorr(residuals(:,i)));
end
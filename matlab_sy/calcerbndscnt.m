%jcsmooth
function [datout]=calcerbndscnt(datin,winsize)
% FFvals % a 1xn vector
nv=length(datin);
% smootherwidth=20; % vary this
% smoothedactual=runave2(FFvals,smootherwidth); % smooth - vary
% this function?
resamplesize=1000;
% clear smoothedtest;
if(length(datin)>2*winsize)
for i=1:resamplesize
   randvec1=sort(ceil(rand(1,nv)*nv));
   datout(i,:)=runmed(datin(randvec1),winsize); % smooth
end
else
    datout=[]
end
% This loop took less than 20 seconds for a 1x141 FFvals vector on my
% computer.
% Results look exactly the same with resamplesize=1000, so reduce resamplesize

% figure
% % figure;plot(smoothedactual,'b')
% hold on;plot(prctile(datout,95),'r')
% hold on;plot(prctile(datout,5),'c')
% hold on;plot(prctile(datout,50),'r')

%takes a running_median of n columns;starts
%at ceil(win/2); ends at n-ceil(win/2);

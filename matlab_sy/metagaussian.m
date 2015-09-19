%metagaussian.m
%this script calls plotgaussian.m

figure
ax=subplot(1,1,1);
xbnds=[-10:.01:10]
mu=1
sigma=1;
col='k'
lstyl='-'
flip=0;

plotgaussian(ax,xbnds,mu,sigma,col,lstyl,1)

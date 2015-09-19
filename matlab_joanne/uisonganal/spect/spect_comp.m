function []= spect_comp(ss1,ss2)


%recover variables

sn1=ss1(1,:);
bf1=ss1(2,:);
mod1=ss1(3,:);
br1=ss1(4,:);
pva1=ss1(5,:);

sn2=ss2(1,:);
bf2=ss2(2,:);
mod2=ss2(3,:);
br2=ss2(4,:);
pva2=ss2(5,:);


%print some stuff to screen



%plot histograms of noisiness


%best ratios
v1=br1;
v2=br2;

%set ceiling size for histogram
nstd=3;  %number of standard deviations above mean to plot
step=.2; %size of step for bins

meanv1=mean(v1);
stdv1=std(v1);
meanv2=mean(v2);
stdv2=std(v2);

ceilv1=meanv1+nstd*stdv1;
ceilv2=meanv1+nstd*stdv2;
mceil=max(ceilv1,ceilv2);
cind1=find(v1>mceil);
v1(cind1)=mceil*ones(size(cind1));
cind2=find(v2>mceil);
v2(cind2)=mceil*ones(size(cind2));


[n1,x1]=hist(v1,[0:step:ceil(max([v1,v2]))]);
[n2,x2]=hist(v2,[0:step:ceil(max([v1,v2]))]);
mn1=n1*100/mean(n1);
mn2=n2*100/mean(n2);


h1=figure;
subplot(2,2,1)
plot(x1,mn1);
hold on
plot(x2,mn2,'m-');
hold off
axis([0,mceil+step,0,max([mn1,mn2])]);
title('best ratios');

%log best ratios
v1=log(br1);
v2=log(br2);

%set ceiling size for histogram
nstd=3;  %number of standard deviations above mean to plot
step=.2; %size of step for bins

meanv1=mean(v1);
stdv1=std(v1);
meanv2=mean(v2);
stdv2=std(v2);

ceilv1=meanv1+nstd*stdv1;
ceilv2=meanv1+nstd*stdv2;
mceil=max(ceilv1,ceilv2);
cind1=find(v1>mceil);
v1(cind1)=mceil*ones(size(cind1));
cind2=find(v2>mceil);
v2(cind2)=mceil*ones(size(cind2));


[n1,x1]=hist(v1,[0:step:ceil(max([v1,v2]))]);
[n2,x2]=hist(v2,[0:step:ceil(max([v1,v2]))]);
mn1=n1*100/mean(n1);
mn2=n2*100/mean(n2);

subplot(2,2,4)
plot(x1,mn1);
hold on
plot(x2,mn2,'m-.');
hold off
legend(num2str(mean(log(br1))),num2str(mean(log(br2))));
axis([0,mceil+step,0,max([mn1,mn2])]);
xlabel('log best ratios');


%peak vs avg
v1=pva1;
v2=pva2;

%set ceiling size for histogram
nstd=3;  %number of standard deviations above mean to plot
step=.2; %size of step for bins

meanv1=mean(v1);
stdv1=std(v1);
meanv2=mean(v2);
stdv2=std(v2);

ceilv1=meanv1+nstd*stdv1;
ceilv2=meanv1+nstd*stdv2;
mceil=max(ceilv1,ceilv2);
cind1=find(v1>mceil);
v1(cind1)=mceil*ones(size(cind1));
cind2=find(v2>mceil);
v2(cind2)=mceil*ones(size(cind2));

[n1,x1]=hist(v1,[0:step:ceil(max([v1,v2]))]);
[n2,x2]=hist(v2,[0:step:ceil(max([v1,v2]))]);
mn1=n1*100/mean(n1);
mn2=n2*100/mean(n2);

subplot(2,2,2)
plot(x1,mn1);
hold on
plot(x2,mn2,'m-.');
hold off
legend(num2str(mean(log(pva1))),num2str(mean(log(pva2))));
axis([0,mceil+step,0,max([mn1,mn2])]);
xlabel('peak vs avg');

%signal to noise
v1=sn1;
v2=sn2;

%set ceiling size for histogram
nstd=3;  %number of standard deviations above mean to plot
step=4; %size of step for bins

meanv1=mean(v1);
stdv1=std(v1);
meanv2=mean(v2);
stdv2=std(v2);

ceilv1=meanv1+nstd*stdv1;
ceilv2=meanv1+nstd*stdv2;
mceil=max(ceilv1,ceilv2);
cind1=find(v1>mceil);
v1(cind1)=mceil*ones(size(cind1));
cind2=find(v2>mceil);
v2(cind2)=mceil*ones(size(cind2));

[n1,x1]=hist(v1,[0:step:ceil(max([v1,v2]))]);
[n2,x2]=hist(v2,[0:step:ceil(max([v1,v2]))]);
mn1=n1*100/mean(n1);
mn2=n2*100/mean(n2);

subplot(2,2,1)
plot(x1,mn1);
hold on
plot(x2,mn2,'m-.');
hold off
axis([0,mceil+step,0,max([mn1,mn2])]);
xlabel('signal to noise (norm?)');



%best freq
v1=bf1;
v2=bf2;

%set ceiling size for histogram
nstd=3;  %number of standard deviations above mean to plot
step=40; %size of step for bins

meanv1=mean(v1);
stdv1=std(v1);
meanv2=mean(v2);
stdv2=std(v2);

ceilv1=meanv1+nstd*stdv1;
ceilv2=meanv1+nstd*stdv2;
mceil=max(ceilv1,ceilv2);
cind1=find(v1>mceil);
v1(cind1)=mceil*ones(size(cind1));
cind2=find(v2>mceil);
v2(cind2)=mceil*ones(size(cind2));

[n1,x1]=hist(v1,[0:step:ceil(max([v1,v2]))]);
[n2,x2]=hist(v2,[0:step:ceil(max([v1,v2]))]);
mn1=n1*100/mean(n1);
mn2=n2*100/mean(n2);

subplot(2,2,3)
plot(x1,mn1);
hold on
plot(x2,mn2,'m-.');
hold off
axis([0,mceil+step,0,max([mn1,mn2])]);
xlabel('best freq');





%figure

h2=figure;
subplot(2,1,1)
plot(sn1,log(br1),'+');
hold
plot(sn2,log(br2),'m+');
hold off
xlabel('signal to noise (before norm)');
ylabel('log(best ratio)');

subplot(2,1,2);
plot(log(br1),bf1,'+');
hold
plot(log(br2),bf2,'+m');
hold off
xlabel('log(best ratio)');
ylabel('best frequency');



delete(h2);
subplot(2,2,1);





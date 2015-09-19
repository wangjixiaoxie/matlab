function  temp_plot(sp1, sp2, range);

%display normal and scrambled spectra over "range" columns  




figure

subplot(2,1,1);
imagesc(flipud(sp1(:,range)));
subplot(2,1,2);
imagesc(flipud(sp2(:,range)));
subplot(2,1,1);

figure
ind=randperm(length(range));
ind=ind+min(range)-1;


subplot(2,1,1);
imagesc(flipud(sp1(:,ind)));
subplot(2,1,2);
imagesc(flipud(sp2(:,ind)));
subplot(2,1,1);


figure
subplot(2,1,1);
m1=mean(sp1);
m2=mean(sp2);
[ymax,i]=max([m1,m2]);
plot(m1(range));
set(gca,'ylim',[0,ymax*1.1]);
subplot(2,1,2);
m2=mean(sp2);
plot(m2(range));
set(gca,'ylim',[0,ymax*1.1]);
subplot(2,1,1);
ylabel('mean amplitude')


%master plot of pitchshiftlesion example.


%first get data together for evren's bird.
subplot(411)

dat=[]
pth{1}='/oriole7/dir1/bk13w63/'
cmd=['cd ' pth{1}];
eval(cmd);
STDN=datenum([2007,5,17]);
BASELINE=2407;
dat(length(dat)+1).fn='trigtest/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='trigtest2/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='ampon/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='ampon2/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='ampon3/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='ampon4/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/screen/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/trigtest/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/trigtest2/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/ampon/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/ampoff/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/ampon2/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ko';
dat(length(dat)+1).fn='postlesion/ampoff2/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='ko';

hold on;
for ii=1:length(dat)
	load(dat(ii).fn);
	sym=dat(ii).sym;
	plt_freqdaymodtw(fv,tbinshft,fbins,BASELINE,sym,STDN,1,[0,25],3,1,3);
end
grid on;hold on;yy=ylim;
title('Bk13W63 Fundamental Freq','FontS',10,'FontW','b');
set(gca,'FontS',10,'FontW','b');
xlabel('Time Since PB Start (days)','FontS',10,'FontW','b');
ylabel('2nd Harmonic Freq (Hz)','FontS',10,'FontW','b');
xx=xlim;xlim(xx+[-1,1]);
saveas(gcf,'FREQ.fig','fig');


%for bk41bk14
ax=subplot(412);
cd /oriole2/bk41bk14/datasum
load pathvals1-analdata.mat
STDN=datenum([2009,5,28])
BASELINE=3432;








plt_freqdaymodtw(avls.fvcomb{1},tbinshft,[3000 4500], BASELINE,sym,STDN,1,[0,25],1,1,1)
% axis([-16 25 3 4])

%CONTROL
ax=subplot(413);
dat=[];
pth{1}='/mnt/removable1/behav_fig_data/pk100bk68/'
cmd=['cd ' pth{1}];
eval(cmd);
STDN=datenum([2006,5,17])+122;
BASELINE=2307;
dat(length(dat)+1).fn='ampoff2/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='bs';
dat(length(dat)+1).fn='ampon3/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ro';
dat(length(dat)+1).fn='ampon3/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='rv';
dat(length(dat)+1).fn='ampoff3/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='bs';
dat(length(dat)+1).fn='ampoff3/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='bs';
dat(length(dat)+1).fn='ampon4/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ro';
dat(length(dat)+1).fn='ampon5/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='r*';
dat(length(dat)+1).fn='ampoff5/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='bs';

hold on;

for ii=1:length(dat)
	load(dat(ii).fn);
	sym=dat(ii).sym;
	plt_freqdaymodtw(fv,tbinshft,fbins,BASELINE,sym,STDN,1,[0,45],3,1,3);
end


ax=subplot(414);
dat=[];
pth{1}='/mnt/removable1/pk30bk79/'
cmd=['cd ' pth{1}];
eval(cmd);
STDN=733400;
BASELINE=2382;
FLIP=-1;
dat(length(dat)+1).fn='prelesion/screen/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='bs';
dat(length(dat)+1).fn='prelesion/ampon/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ro';
dat(length(dat)+1).fn='prelesion/ampon2/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='rd';
dat(length(dat)+1).fn='prelesion/ampon3/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='ro';
dat(length(dat)+1).fn='prelesion/ampon4/BATCHTRAINCATCH_FFVN4';
dat(length(dat)).sym='rd';
dat(length(dat)+1).fn='prelesion/ampoff/BATCHTRAIN_FFVN4';
dat(length(dat)).sym='bs';


hold on;

for ii=1:length(dat)
	load(dat(ii).fn);
	sym=dat(ii).sym;
	plt_freqdaymodtw(fv,tbinshft,fbins,BASELINE,sym,STDN,1,[-15,45],3,1,3,FLIP);
end
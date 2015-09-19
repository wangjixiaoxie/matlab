%safplotscript
%used for plotting results of saffeedback expts.
figure
fn='bk68bk82_170507_1552.3.cbin';

[dat,fs]=ReadCbinFile(fn);
rd=readrecf(fn);
tt=rd.ttimes;
sdt=dat(:,1);
pbdt=dat(:,2);
fs=32000
figure
ax4(1)=subplot(2,1,1)


evspect(sdt,fs);
hold on;
xvla=[tt/1000 tt/1000+.06];
yvla=[500 500];
plot(xvla,yvla,'k', 'Linewidth', 5)

ax4(2)=subplot(2,1,2)
sdt=sdt-mean(sdt);
wndat=sdt;
for ii=1:length(tt)
    
    wndat(floor(tt(ii)/1000*fs):floor(tt(ii)/1000*fs+.06*fs)-1)=wndat(floor(tt(ii)/1000*fs):floor(tt(ii)/1000*fs+.06*fs)-1)+100000*datfile;
end
evspect(wndat,fs);
linkaxes(ax4)
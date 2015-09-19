function evautolbl_screen(od,CHANSPEC);
% evautolbl_screen(od,CHANSPEC);

if (~exist('CHANSPEC'))
	CHANSPEC='obs0';
end

fig=figure;
set(fig,'Position',[5,508,1267,420]);

for ii=1:length(od)
	[dat,fs]=evsoundin('',od(ii).fn,CHANSPEC);
	[sm,sp,t,f]=evsmooth(dat,fs,100);
	ax1=subplot(211);imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
	tt=0.5*(od(ii).onsets+od(ii).offsets)*1e-3;
	ax2=subplot(212);plot(tt,od(ii).mvals,'bs-');grid on;
	linkaxes([ax1,ax2],'x');
	xlim([0,t(end)]);
	pan xon;zoom on;
	drawnow;
	pause;
end
return;

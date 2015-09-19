%needs odat ii CS
[dat,fs]=evsoundin('',odat(ii).fn,CS);
[sm,sp,t,f]=evsmooth(dat,fs,100);

ax1=subplot(211);imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
x=0.5e-3*(odat(ii).onsets+odat(ii).offsets);
ax2=subplot(212);plot(x,odat(ii).mvals,'bs-')
grid on;
linkaxes([ax1,ax2],'x');
xlim([0,t(end)])
pan xon;zoom xon;

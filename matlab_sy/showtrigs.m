function showtrigs(batchfile,SPTH);
%

if (~exist('SPTH'))
	SPTH=0.01;
end

fid=fopen(batchfile,'r');
while (1)
	fn = fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	[sng,fs]=readevtaf(fn,'0r');
	[trg,fs]=readevtaf(fn,'1r');
	[sm,sp,t,f]=evsmooth(sng,fs,SPTH);

	imagesc(t,f,log(abs(sp)));imstyle;hold on;grid on;
	trg = trg*1e4/max(trg);
	plot([1:length(trg)]/fs,trg,'r-');
	title(normstr(fn));
	drawnow;pause;hold off;
end
fclose(fid);
return;

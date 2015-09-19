%templ=load('r93bu81_templates3.dat'); 
[dat,fs]=readevtaf(fn,'0r');
[sm,sp,t,f]=evsmooth(dat,fs,0.02);
nfft=size(templ,1)*2;

[out,ff,tt]=specgram(dat,nfft,fs,hanning(nfft),0);
out = abs(out);
out=out(1:nfft/2,:);

df = zeros([length(out),size(templ,2)]);
for ii = 1:length(out)
	%out(:,ii) = out(:,ii) - min(out(:,ii));
	%out(:,ii) = out(:,ii)./max(out(:,ii));
    outn = out(:,ii);outn(1:6)=0.0;
	outn = outn./sqrt(outn.'*outn);

	for jj = 1:size(templ,2)
		templn = templ(:,jj)./sqrt(templ(:,jj).'*templ(:,jj));
		%df(ii,jj) = sum((outn - templn).^2);
		df(ii,jj) = acos(outn.'*templn);
	end
end

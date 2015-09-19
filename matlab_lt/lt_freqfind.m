% 08/23/13 - LT edited so that load on line 27 is just fn, not fn with
% .not.mat added on end.

function vals=freqfind(bt,freqrange,NT,NFFT,USEX);
% vals=freqfind(bt,freqrange,NT,NFFT);
% bt - batch file
% freqrange - [min,max] freq in Hz
% NT = the note you are looking for
% NFFT - number of data points in the template
% vals - output frequency values 
%USEX - if O use regular x file.
% example : vals=frqfind('batch',[7000,8000],'a',128);

if (~exist('NFFT'))
    NFFT=256;
else
    NFFT=NFFT*2;
end
if(~exist('USEX'))
    USEX=1
end

vals=[];
ff=load_batchf(bt);
for ii=1:length(ff)
	fn=ff(ii).name;
	rd=readrecf(fn,USEX);
	load(fn);
    disp(fn);
    
	[dat,fs]=evsoundin('',fn,'obs0');

	if (~isfield(rd,'ttimes'))
		tt=[];
	else
		tt=rd.ttimes*1e-3;
	end
	for jj=1:length(tt)
		pp=find((onsets<=tt(jj)*1e3)&(offsets>=tt(jj)*1e3));
		if (~strcmp(labels(pp),NT))
			continue;
		end
		inds=fix(tt(jj)*fs)+[-(NFFT-1):0];
		dat2=dat(inds);
		fdat=abs(fft(hamming(NFFT).*dat2));

		ffv=get_fft_freqs(NFFT,32e3);
		ffv=ffv(1:end/2);
		
		tempvals=[];
        fbins=freqrange(1,:);
		for kk=1:size(fbins,1)
			inds2=find((ffv>=fbins(kk,1))&(ffv<=fbins(kk,2)));
			[y,i]=max(fdat(inds2));
			i=i+inds2(1)-1;
			i=i+[-3:3];
			tempvals=[tempvals,sum(ffv(i).*fdat(i).')./sum(fdat(i))];
		end
		vals=[vals;tempvals];
	end
end

return;
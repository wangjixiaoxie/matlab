function vals=evtaf_freq(bt,fbins,NT,NFFT,CS,USEX);
% vals=evtaf_freq(bt,fbins,NT,NFFT,CS,USEX);
% vals=[vals;fn2datenum(fn),tempvals,IND,filenum];
%
%  bt - batch file
%  fbins - [Min Freq, Max Freq] to search for peaks (in Hz)
%  NT - target note
%  NFFT - length of template 
%  CS - chan spec
%  USEX - if == 1 look in X.rec for trigger times
% returns vals - 
%   vals=[datenum of file , FREQ vals , Note Index , file number];
% usage example:
%  vals=evtaf_freq('batch.train',[5000,6000],'a',128,'obs0',0);

Nbins=3;
NFFT=NFFT*2;
vals=[];
ff=load_batchf(bt);
for ii=1:length(ff)
	fn=ff(ii).name;
	ppp=findstr(fn,'.cbin');
	pppp=findstr(fn,'.');
	tmp=find(pppp<ppp);pppp=pppp(tmp(end));
	filenum=str2num(fn(pppp+1:ppp-1));
	if (~exist(fn,'file'))
		continue;
	end
	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	rd=readrecf(fn,USEX);
	load([fn,'.not.mat']);

	[dat,fs]=evsoundin('',fn,CS);

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
		IND=pp(1);
		inds=fix(tt(jj)*fs)+[-(NFFT-1):0];
		dat2=dat(inds);
		fdat=abs(fft(hamming(NFFT).*dat2));

		ffv=get_fft_freqs(NFFT,fs);
		ffv=ffv(1:end/2);
		
		tempvals=[];
		for kk=1:1%size(fbins,1)
			inds2=find((ffv>=fbins(kk,1))&(ffv<=fbins(kk,2)));
			[y,i]=max(fdat(inds2));
			i=i+inds2(1)-1;
			i=i+[-Nbins:Nbins];
			tempvals=[tempvals,sum(ffv(i).*fdat(i).')./sum(fdat(i))];
		end
		vals=[vals;fn2datenum(fn),tempvals,IND,filenum];
	end
end
return;

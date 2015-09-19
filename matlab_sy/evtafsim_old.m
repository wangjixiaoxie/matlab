function [cnt1,cnt2,FullSpec,th] = ...
           evtafsim(rawsong,rs_Fs,thresh,RefSpec1,RefSpec2,ASpecFile);
% returns the threshold in case it was set in the program
% simtaf with the clocks and multiple templates put in
%

% have no idea what this number is
% seems like the normilaztion would make this not be required
MultFac = .92;

zbins= 6; % this is what TAF does i guess
fs   = 22050;
blen = 128;
nfft = blen*2;

%[rawsong,rs_Fs]=wavread(birdsongfile);
rawsong = resample(rawsong,fs,rs_Fs);

%if (exist(RefF1,'file'))
%	RefSpec1 = load(RefF1);
%else
%	disp(['Ref File 1 does not exists : ',RefF1]);
%end

%if (exist(RefF2,'file'))
%	RefSpec2 = load(RefF2);
%else
%	disp(['Ref File 2 does not exists : ',RefF2]);
%end

if (exist('ASpecFile'))
	if (exist(ASpecFile,'file'))
		AmSp = load(ASpecFile);
	else
		AmSp = 0.0;
	end
else
	AmSp = 0.0;
end


rsongu = MultFac.*(rawsong - mean(rawsong));
hammy  = hamming(2*blen);
kern   = boxkern(3,3); % looks like the convolution of 2 box kernels

nrep = floor(length(rsongu)/blen)-1;
vals1 = zeros([nrep,1]);vals2 = zeros([nrep,1]);
FullSpec = zeros([nrep,blen]);

for ii = 1:nrep
	ind1 = (ii-1)*blen + 1;
	ind2 = ind1 + 2*blen - 1;
	datchunk = rsongu(ind1:ind2);
	fdatchunk = abs(fft(hammy.*datchunk));
	fdatchunk = fdatchunk(1:blen) - AmSp;
	fdatchunk = smooth(fdatchunk,kern);
	fdatchunk(1:zbins) = 0.0;
	sp = fdatchunk;
	sp(find(sp)<0) = 0.0;
	sp = sp - min(sp);
	sp = sp./max(sp);

	sp_comp = (sp - RefSpec1);
	vals1(ii) =  sp_comp.'*sp_comp;

	sp_comp = (sp - RefSpec2);
	vals2(ii) =  sp_comp.'*sp_comp;

	FullSpec(ii,:) = sp.';
end

if (length(thresh)==1)
	thresh = [1;1]*thresh;
end
if (length(thresh)==0)
	figure;
	subplot(211);imagesc(FullSpec.');
	set(gca,'YD','no');m=colormap('gray');colormap(m(end:-1:1,:));
	subplot(212);plot(vals1,'b');hold;grid;plot(vals2,'r');
	drawnow;
	disp(['Hit return to get Threshold set cross hairs!']);pause;
	thresh = zeros([2,1]);
	disp(['Pick threshold 1 - BLUE LINE']);
	g=ginput(1);thresh(1)=g(2);
	disp(['Pick threshold 2 - RED  LINE']);
	g=ginput(1);thresh(2)=g(2);
end
th = thresh;
%disp(['Using threshold values of :',num2str(th(1)),' and ',num2str(th(2))]);


cnt1 = zeros([nrep,1]);cnt2=cnt1;
cnt1(1:zbins) = [1:zbins].';cnt2(1:zbins) = [1:zbins].';
for ii = (zbins+1):nrep
	if (vals1(ii)<=thresh(1))
		cnt1(ii) = 0;
	else
		cnt1(ii) = cnt1(ii-1) + 1;
	end

	if (vals2(ii)<=thresh(2))
		cnt2(ii) = 0;
	else
		cnt2(ii) = cnt2(ii-1) + 1;
	end
end

return;

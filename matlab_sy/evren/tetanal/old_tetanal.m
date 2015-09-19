%[dat,fs]=evsoundin('',fn,'obs1');
%[tdat,fs]=evsoundin('',fn,'obs2');
%dat=[dat,tdat];
%[tdat,fs]=evsoundin('',fn,'obs3');
%dat=[dat,tdat];
%[tdat,fs]=evsoundin('',fn,'obs4');
%dat=[dat,tdat];



[data,fs]=ReadCbinFile(fname);


%TH=-2000;
refrac=0.25e-3;
NChan=4;

pk_win_sz=ceil(refrac*fs);
for iChan = 1:NChan
	spki=find((dat(1:end-1,iChan)>TH)&(dat(2:end,iChan)<=TH))+1;
	spkt=spki/fs;
	disp(length(spkt));

	indsv=zeros([length(spkt),1]);
	ind1=1;ind2=2;
	nindsv=0;
	while (ind1<length(spkt))
		if (spkt(ind2)-spkt(ind1)<refrac)
			ind2=ind2+1;
		else
			nindsv=nindsv+1;
			indsv(nindsv)=ind1;
			ind1=ind2;
			ind2=ind1+1;
		end
		if (ind2>length(spkt))
			break;
		end
	end
	spki = fix(round(spkt(indsv(1:nindsv))*fs));

	% take min over the next 0.5 ms as the peak
	spkt=0*spki;spkamp=0*spki;
	for ii=1:length(spki)
		if (spki(ii)<=length(dat)-pk_win_sz)
			[vmn,imn]=min(dat(spki(ii)+[0:pk_win_sz],iChan));
		else
			[vmn,imn]=min(dat(spki(ii):end,iChan));
		end
		spkt(ii)=(imn-1+spki(ii))/fs;
		%spkamp(ii)=abs(vmn);
	end
	spki = fix(round(spkt*fs));

	disp(length(spkt));
	spkt_sv{iChan}=spkt;

	spkamp = zeros([size(spki,1),NChan]);
	for ii=1:length(spki)
		spkamp(ii,:) = -1*(min(dat(spki(ii)+[-2:2],:)));
	end
	spka_sv{iChan}=spkamp;
end

function [adat1,adat2]=get_ampcont(ff,NT,CS,tbs1,nfft1,dt_note);
% [adat1,adat2]=get_ampcont(ff,NT,CS,tbs1,nfft1,dt_note);
%

adat1=[];
for kk=1:length(tbs1)
        adat1(kk).cnt=0;
        adat1(kk).vals=[];
end

adat2=[];
for kk=1:size(dt_note,1)
        %control
        adat2(kk).cnt=0;
        adat2(kk).mnt=0;
        adat2(kk).stt=0;
        adat2(kk).vals=[];
        adat2(kk).badcnt=0;
end

FREQRNG=[100,1e4];

for ii=1:length(ff)
	[dat,fs]=evsoundin('',ff(ii).name,CS);
	load([ff(ii).name,'.not.mat']);
	sm=evsmooth(dat,fs,0);
	%get norm
	tmpinds=[floor(onsets(1)*1e-3*fs),ceil(offsets(end)*1e-3*fs)];
	normval=mean(sqrt(sm(tmpinds(1):tmpinds(2))));
	DT=1./fs;

	pl=findstr(labels,NT);
	for jj=1:length(pl)
		%data around the target note
		ons=onsets(pl(jj))*1e-3;
		for kk=1:length(tbs1)
			inds=floor((ons+tbs1(kk))*fs)+[0:nfft1(kk)-1];
			dattmp=sqrt(sm(inds));

			tmpval=sum(dattmp).*DT./normval;
			adat1(kk).vals=[adat1(kk).vals;tmpval];
			adat1(kk).cnt=adat1(kk).cnt+1;
		end

		%neighboring notes
		for kk=1:size(dt_note,1)
			dts=dt_note(kk,1:2);
			tbs2=dt_note(kk,3);
			nfft2=fix(dt_note(kk,4));
			pp=find((onsets*1e-3>=ons+dts(1))&...
				(onsets*1e-3<=ons+dts(2)));
			if (length(pp)==0)
				adat2(kk).badcnt=adat2(kk).badcnt+1;
				continue;
			end

			ons2=onsets(pp(1))*1e-3;
			inds=floor((ons2+tbs2)*fs)+[0:nfft2-1];
			dattmp=sqrt(sm(inds));

			adat2(kk).cnt=adat2(kk).cnt+1;
			adat2(kk).mnt=adat2(kk).mnt+(ons2-ons);
			adat2(kk).stt=adat2(kk).stt+(ons2-ons).^2;

                        tmpval=sum(dattmp).*DT./normval;
                        adat2(kk).vals=[adat2(kk).vals;tmpval];
		end
	end
end
for kk=1:length(adat2)
	adat2(kk).mnt=adat2(kk).mnt./adat2(kk).cnt;
	adat2(kk).stt=sqrt((adat2(kk).stt./adat2(kk).cnt)-(adat2(kk).mnt.^2));
end
return;

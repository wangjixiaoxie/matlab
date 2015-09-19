function [adat]=get_ampcont2(ff,NT,CS,dt_note,TH,MININT,MINDUR);
% [adat]=get_ampcont2(ff,NT,CS,dt_note,TH,MININT,MINDUR);
%

adat=[];
for kk=1:size(dt_note,1)
        %control
        adat(kk).cnt=0;
        adat(kk).mnt=0;
        adat(kk).stt=0;
        adat(kk).vals=[];
        adat(kk).badcnt=0;
end

for ii=1:length(ff)
	[dat,fs]=evsoundin('',ff(ii).name,CS);
	load([ff(ii).name,'.not.mat']);
	sm=evsmooth(dat,fs,0);
	sm=sqrt(sm);
	%get norm
	tmpinds=[floor(onsets(1)*1e-3*fs),ceil(offsets(end)*1e-3*fs)];
	normval=mean(sm(tmpinds(1):tmpinds(2)));
	DT=1./fs;

	smn=sm./normval;
	[nons,noffs]=evsegment(smn,fs,MININT,MINDUR,TH);

	pl=findstr(labels,NT);
	for jj=1:length(pl)
		[y,i]=min(abs(onsets(pl(jj))-nons));
                ons=nons(i(1))*1e-3;

		%neighboring notes
		for kk=1:size(dt_note,1)
			dts=dt_note(kk,1:2);
			pp=find((nons*1e-3>=ons+dts(1))&...
				(nons*1e-3<=ons+dts(2)));
			if (length(pp)==0)
				adat(kk).badcnt=adat(kk).badcnt+1;
				continue;
			end

			ons2=nons(pp(1))*1e-3;
			inds=floor(ons2*fs)+[0:fix(dt_note(kk,3)*fs)];
			dattmp=sm(inds);

			adat(kk).cnt=adat(kk).cnt+1;
			adat(kk).mnt=adat(kk).mnt+(ons2-ons);
			adat(kk).stt=adat(kk).stt+(ons2-ons).^2;

                        tmpval=sum(dattmp).*DT./normval;
                        adat(kk).vals=[adat(kk).vals;tmpval];
		end
	end
end
for kk=1:length(adat)
	adat(kk).mnt=adat(kk).mnt./adat(kk).cnt;
	adat(kk).stt=sqrt((adat(kk).stt./adat(kk).cnt)-(adat(kk).mnt.^2));
end
return;

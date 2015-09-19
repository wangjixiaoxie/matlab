function [adat1,adat2]=get_ffcont(ff,NT,CS,tbs1,nfft1,dt_note,fbins1,fbins2);
% [adat1,adat2]=get_ffcont(ff,NT,CS,tbs1,nfft1,dt_note);
%

adat1=[];
for kk=1:length(tbs1)
        adat1(kk).avfdat=zeros([nfft1(kk)/2,1]);
        adat1(kk).stfdat=zeros([nfft1(kk)/2,1]);
        adat1(kk).cnt=0;
        adat1(kk).vals=[];
end

adat2=[];
for kk=1:size(dt_note,1)
        nfft2=fix(dt_note(kk,4));
        %control
        adat2(kk).avfdat=zeros([nfft2/2,1]);
        adat2(kk).stfdat=zeros([nfft2/2,1]);
        adat2(kk).cnt=0;
        adat2(kk).mnt=0;
        adat2(kk).stt=0;
        adat2(kk).vals=[];
        adat2(kk).badcnt=0;
end

NPNTS=2;disp(['NPNTS = ',num2str(NPNTS)]);

for ii=1:length(ff)
	[dat,fs]=evsoundin('',ff(ii).name,CS);
	load([ff(ii).name,'.not.mat']);
	pl=findstr(labels,NT);
	for jj=1:length(pl)
		%data around the target note
		ons=onsets(pl(jj))*1e-3;
		for kk=1:length(tbs1)
			inds=floor((ons+tbs1(kk))*fs)+[0:nfft1(kk)-1];
			dattmp=dat(inds);
			fdatmp=abs(fft(dattmp.*hamming(length(dattmp))));
			fdatmp=fdatmp(1:end/2);
			adat1(kk).avfdat=adat1(kk).avfdat+fdatmp;
			adat1(kk).stfdat=adat1(kk).stfdat+fdatmp.^2;
			adat1(kk).cnt=adat1(kk).cnt+1;

			ffv=get_fft_freqs(length(fdatmp)*2,fs);
			ffv=ffv(1:end/2);
			fbins=fbins1(kk).fb;
			mxvals=zeros([size(fbins,1),1]);
			for mm=1:size(fbins,1)
				pp=find((ffv>=fbins(mm,1))&(ffv<=fbins(mm,2)));
				[yyy,mxi]=max(fdatmp(pp));
				mxi=mxi+pp(1)-1;
				mxi=mxi+[-NPNTS:NPNTS];
				mxvals(mm)=sum(ffv(mxi).'.*fdatmp(mxi)./sum(fdatmp(mxi)));
			end
			tmpval=mean(mxvals.'./[1:length(mxvals)]);
			adat1(kk).vals=[adat1(kk).vals;tmpval];
		end

		%neighboring notes
		for kk=1:size(dt_note,1)
			dts=dt_note(kk,1:2);tbs2=dt_note(kk,3);nfft2=fix(dt_note(kk,4));
			pp=find((onsets*1e-3>=ons+dts(1))&(onsets*1e-3<=ons+dts(2)));
			if (length(pp)==0)
				adat2(kk).badcnt=adat2(kk).badcnt+1;
				continue;
			end
			ons2=onsets(pp(1))*1e-3;
			inds=floor((ons2+tbs2)*fs)+[0:nfft2-1];
			dattmp=dat(inds);
			fdatmp=abs(fft(dattmp.*hamming(length(dattmp))));
			fdatmp=fdatmp(1:end/2);

			adat2(kk).avfdat=adat2(kk).avfdat+fdatmp;
			adat2(kk).stfdat=adat2(kk).stfdat+fdatmp.^2;
			adat2(kk).cnt=adat2(kk).cnt+1;
			adat2(kk).mnt=adat2(kk).mnt+(ons2-ons);
			adat2(kk).stt=adat2(kk).stt+(ons2-ons).^2;

                        ffv=get_fft_freqs(length(fdatmp)*2,fs);
                        ffv=ffv(1:end/2);
                        fbins=fbins2(kk).fb;
                        mxvals=zeros([size(fbins,1),1]);
                        for mm=1:size(fbins,1)
                                pp=find((ffv>=fbins(mm,1))&(ffv<=fbins(mm,2)));
                                [yyy,mxi]=max(fdatmp(pp));
                                mxi=mxi+pp(1)-1;
                                mxi=mxi+[-NPNTS:NPNTS];
                                mxvals(mm)=sum(ffv(mxi).'.*fdatmp(mxi)./sum(fdatmp(mxi)));
                        end
                        tmpval=mean(mxvals.'./[1:length(mxvals)]);
                        adat2(kk).vals=[adat2(kk).vals;tmpval];
		end
	end
end
for kk=1:length(adat1)
	adat1(kk).avfdat=adat1(kk).avfdat./adat1(kk).cnt;
	adat1(kk).stfdat=adat1(kk).stfdat./adat1(kk).cnt...
		       - adat1(kk).avfdat.^2;
	adat1(kk).stfdat=sqrt(adat1(kk).stfdat);
end
for kk=1:length(adat2)
	adat2(kk).avfdat=adat2(kk).avfdat./adat2(kk).cnt;
	adat2(kk).stfdat=adat2(kk).stfdat./adat2(kk).cnt...
		       - adat2(kk).avfdat.^2;
	adat2(kk).stfdat=sqrt(adat2(kk).stfdat);
	
	adat2(kk).mnt=adat2(kk).mnt./adat2(kk).cnt;
	adat2(kk).stt=sqrt((adat2(kk).stt./adat2(kk).cnt)-(adat2(kk).mnt.^2));
end
return;

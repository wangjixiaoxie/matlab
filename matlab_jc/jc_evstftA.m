function FFOUT=jc_evstftA(shiftedCOMPARISON,mini,maxi,windowsize,stepsize)



fs=32e3;




nframe=windowsize;
overlap=windowsize/stepsize;
overfac=1./overlap;
stepsize=nframe/overlap;

% these are the bin frequencies
%ff=get_fft_freqs(nframe,fs);
%ff=ff(1:((nframe/2)+1)).';

%expected phase advance
freqperbin=fs/nframe;
%dphi_exp = 2*pi*ff.*(nframe/overlap)./fs;
dphi_exp = 2*pi*[0:((nframe/2)-1)].'*stepsize/nframe;
wind=hanning(nframe);


FFOUT=[];
for lll=1:size(shiftedCOMPARISON,1)
	st_ind=1;en_ind=st_ind+nframe-1;
	phvalsv=zeros([(nframe/2),1]);%start all phases at zero
	ffvsv=[];phsv=[];dphsv=[];powsv=[];
	spgrm=[];tvals=[];ffdev_sv=[];
	cnt=1;
	dat2=shiftedCOMPARISON(lll,:).';
	while (1)
		if (en_ind>length(dat2))
			break;
		end
		cnt=cnt+1;
	
		dat_tmp=wind.*dat2(st_ind:en_ind);
		fdat=fft(dat_tmp);
		fdat=fdat(1:(nframe/2));
	
		ph=atan2(imag(fdat),real(fdat));
		amp=(real(fdat).^2 + imag(fdat).^2);

		delt_ph = (ph-phvalsv) - dphi_exp;
		%delt_ph = delt_ph - 2*pi*round(delt_ph./pi./2);
		%below is from SB code
		NPI = fix(delt_ph/pi);
		if (NPI >= 0) 
			NPI = NPI + mod(NPI,2);
		else
			NPI = NPI - mod(NPI,2);
		end
		delt_ph =delt_ph - pi*NPI;
	
		ffdev=delt_ph*overlap/2/pi;
		actff=(fs/nframe)*(([0:(nframe/2)-1].')+ffdev);

		ffdev_sv=[ffdev_sv,ffdev];
	
		ffvsv = [ffvsv,actff];
		phsv  = [phsv,ph];
		dphsv = [dphsv,delt_ph];
		powsv = [powsv,amp];

		phvalsv=ph;

		st_ind=st_ind+stepsize;
		en_ind=st_ind+nframe-1;
		tvals=[tvals;0.5*(en_ind+st_ind)/fs];
	end

	ffsv=[];
	fbins=[mini,maxi;mini*2,maxi*2;mini*3,maxi*3];
	disp('Hey FBINS');
	for ii=1:size(ffvsv,2)
		tmp=ffvsv(:,ii);
		tmp2=[];
		%for jj=1:size(fbins,1)
		for jj=1:size(fbins,1)
                ppp=find((tmp>fbins(jj,1))&(tmp<=fbins(jj,2)));
                if (length(ppp)>0)
                    [y,i]=max(powsv(ppp,ii));
                    inds=ppp(i);%+[-1:1];
                    inds=inds(find((inds>0)&(inds<=length(ffvsv(:,ii)))));
                    tmp2=[tmp2,sum(ffvsv(inds,ii).*powsv(inds,ii))./sum(powsv(inds,ii))./jj];
                    %tmp2=ffvsv(inds,ii);
                end
                
                
                
                
		end
		ffsv=[ffsv;mean(tmp2)];
	end
	
	FFOUT=[FFOUT;ffsv.'];
end

%figure(1);hold on;plot(tvals,FFOUT(1,:),'r.-');
for ii=1:size(shiftedCOMPARISON,1)
	[sm,sp,t,f]=evsmooth(shiftedCOMPARISON(ii,:),32e3,0.01);
	%figure(11);cla;   
	%imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
	%hold on;grid on;plot(tvals,FFOUT(ii,:),'r-');
	%cmgray;
	%drawnow;pause;
end
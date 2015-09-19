function [trigs,ff2val]=jctafdelay2(batchtmp)
% does simulate evtaf_freq_mult_delay
% capability of using multiple templates/counters to recognize note, check
% FFT at the recognition point and then check the FFT n points in the
% future and trigger white noise at the second FFT check point

% eventually, make many of these parameters into input variables
ff=load_batchf(batchtmp);
fs=32000;
NFFT=256;
cntrng(1).MIN=2;
cntrng(2).MIN=5;
distance=cntrng(2).MIN-cntrng(1).MIN;
cntrng(1).TH=2.8;
cntrng(2).TH=2.8;
conting=[2436 2416];
hitbelow=[1 0];
Nbins=3;
mmm=1;
fbins=[2200 2800];
refracperiod=100; % in milliseconds
refracwidth=refracperiod/8;
for ll=1:length(ff)
    count=0;
    tmpfile=load([ff(ll).name '.tmp']);
    % split tmpfile into a separate vector for each template
        indices=[0:1:-1+length(tmpfile)/length(cntrng)];
        tmpvals=[];
        for ii=1:length(cntrng)
            tmpvals(ii,:)=tmpfile(ii+indices*length(cntrng));
        end
    soundfile=evsoundin('',[ff(ll).name '.cbin'],'obs0');
    soundfile=soundfile(fs*2+1:length(soundfile));
    %%%% Calculate the counters
    counter=zeros(length(cntrng));
    for jjj=1:size(tmpvals,2)
        for kk=1:size(tmpvals,1)
            if tmpvals(kk,jjj)<cntrng(kk).TH
                counter(kk)=counter(kk)+1;
            else
                counter(kk)=0;
            end
            reccounter(kk,jjj)=counter(kk);
        end
    end
    refraconset=0-refracwidth;
    for j=1:size(tmpvals,2)
        % If wn isn't playing
        if j>refraconset+refracwidth
        % Does the note match?
        if reccounter(1,j)==cntrng(1).MIN
            % Does the front of the note match the template?
            % If so, check the first frequency contingency
                dat1=soundfile((j-1)*(NFFT)+1:(j)*(NFFT)); %%%%%
                fdat1=abs(fft(hamming(NFFT).*dat1));
                ffv=get_fft_freqs(NFFT,fs);
                ffv=ffv(1:end/2);
                tempvals1=[];
                for kk=1:1%size(fbins,1)
                    inds2=find((ffv>=fbins(kk,1))&(ffv<=fbins(kk,2)));
                    [y1,i1]=max(fdat1(inds2));
                    i1=i1+inds2(1)-1;
                    i1=i1+[-Nbins:Nbins];
                    tempvals1=[tempvals1,sum(ffv(i1).*fdat1(i1).')./sum(fdat1(i1))];
                end
            % If the first frequency contingency is met, check at the second
            % point in the note.
            if (hitbelow(1)==1 && tempvals1<conting(1)) || (hitbelow(1)==0 && tempvals1>conting(1))
                % Does it match the template?
                if 1 %reccounter(2,j+distance)==cntrng(2).MIN
                    %%%%%%%%%%%%% THIS IS WHERE THE PROBLEM IS
                    %%%%%%%%%%%%% %%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % If so, check the frequency
                    dat1=soundfile((j-1+distance)*(NFFT)+1:(j+distance)*NFFT);
                    fdat1=abs(fft(hamming(NFFT).*dat1));
                    ffv=get_fft_freqs(NFFT,fs);
                    ffv=ffv(1:end/2);
                    tempvals2=[];
                    for kk=1:1%size(fbins,1)
                        inds2=find((ffv>=fbins(kk,1))&(ffv<=fbins(kk,2)));
                        [y1,i1]=max(fdat1(inds2));
                        i1=i1+inds2(1)-1;
                        i1=i1+[-Nbins:Nbins];
                        tempvals2=[tempvals2,sum(ffv(i1).*fdat1(i1).')./sum(fdat1(i1))];
                        ff2val(mmm)=tempvals2;
                        mmm=mmm+1;
                    end
                    if (hitbelow(2)==1 && tempvals2<conting(2)) || (hitbelow(2)==0 && tempvals2>conting(2))
                        count=count+1;
                        trigs(ll).times(count)=2+((NFFT)/fs)*(j+distance);
                        refraconset=j+distance;
                    end
                end
            end
        end
        end
        end
    end
end
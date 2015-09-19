function featvect=feature_vect(syllable,Fs)

% featvect=FEATURE_VECT(syllable)
%
% Returns a row vector which is the feature vector of a particular syllable.
%

   %Fs=32000;

   Dur=length(syllable)/Fs;  % Duration of the syllable
   if ~(length(syllable)>1)
       disp('ERROR: EMPTY SYLLABLE IN FEATURE_VECT!!')
       return
   end
   
   % Calculate the loudness of the syllable (the squared amplitude, smoothed with a 2 ms time window)
   window=ceil(0.002*Fs);
   kern=[];
   for i=1:window;  %  A broad gaussian window
       kern(i)=exp(-(i-window/2)*(i-window/2)/(2.5*window));
   end   
   loud=smooth(syllable.*syllable,kern')/sum(kern);  % Calculate loudness

   % Find the peak loudness and record it's amplitude and timing
   [peak,peakat]=max(loud);
   peakat=peakat/length(loud);
   
   % Normalize;
   loud=loud/sum(loud);
   
   % Calculate the spectral density of the syllable
   [sd,sdf]=specgram(syllable,length(syllable),Fs);
   sd=sd(find(sdf<8000));    % Above 8000 Hz, filtered out anyway
   sdf=sdf(find(sdf<8000));
   sd=abs(sd);  % Throw away relative phase information
   % Normalize
   sd=sd/sum(sd);
   
   % Calculate the spectral density of each half of the syllable
   [sds,sdfs]=specgram(syllable,1+ceil(length(syllable)/2),Fs);
	       sd1=abs(sds(:,1));
	       sd2=abs(sds(:,2));
	       sum1=sum(sd1);
	       sum2=sum(sd2);
	       f1=sum(sd1.*sdfs)/sum1;
	       f2=sum(sd2.*sdfs)/sum2;

   % Calculate the peak and spectral entropy of the full spectrogram
   [spec,specf]=specgram(syllable,256,Fs,205);
   for i=1:size(spec,1)*size(spec,2)
	unwrap(i)=abs(spec(i));
   end
   [y,index]=sort(unwrap);
   spec_peak=mean(unwrap(index>0.99*length(unwrap)));
   unwrap=unwrap/sum(unwrap);
   spec_ent=-sum(unwrap.*log2(unwrap))/log2(length(unwrap));

   % Calculate the time at which half of the loudness has occurred
	       cumloud=cumsum(loud);
	       tmp=find(cumloud>0.5);
	       halftime=tmp(1)/length(loud);

   % Calculate the max and circular average of the fractional part of log2(f),
   % weighted appropriately.   -- the fractional part of log2(f) will be
   % the same for the first, second and fourth harmonics.
   circa=mod(log2(sdf(2:end)),1);
   circf=exp(2i*pi*circa);
   circp=sd(2:end)./sdf(2:end);
   circ_avg=0.5+angle(sum(circf.*circp))/(2*pi);
   circ_mag=abs(sum(circf.*circp));
   chis=zeros(30,1);
   for i=1:length(circf)
	 tmp=ceil(30*circa(i));
         if tmp==0
	      tmp=30;
         end
	 chis(tmp)=chis(tmp)+circp(i);
   end
   [tmp,cmax]=max(chis);
   circ_max=cmax/30;

   % Calculate the Cepstrum Peak
	       [cp,cpf]=specgram(sd,length(sd));
	       cp=abs(cp);
	       cp=cp/sum(cp);
	       kern=[];
	       for i=1:20
	            kern(i)=10-abs(10-i);
               end
	       spt=find(smooth(diff(cp),kern')>0);
               if ~isempty(spt)
                   [tmp,cep_peakat]=max(cp(spt(1):end));
                   cep_peakat=cep_peakat+spt(1)-1;
               else
                   cep_peakat=0;
               end

   % Output the Relevant Features
   featvect(1)=sum(sd.*sdf)/100;    % First: The Mean Frequency
   featvect(2)=-100*sum(sd.*log2(sd))/log2(length(sd)); % Second: Entropy of the spectral density (range: 0 to 100; 100 = white noise; 0 = pure note;
   %featvect(3)=10*log10(peak);                  % Third: Peak Loudness
   featvect(3)=10*log2(Dur*1000);             % Fourth: Duration of Syllable (in ms);
   featvect(4)=-100*sum(loud.*log2(loud))/log2(length(loud));    % Fifth: Entropy of the loudness vs. time (range: 0 to 100; 100 = flat amplitude; 0 = one very sharp peak)
   %featvect(6)=100*circ_avg;     %  A harmonic-corrected mean frequency (2^n+s^fv=fundamental frequency)
   %featvect(7)=100*circ_max;     %  A harmonic-corrected max frequency
   featvect(5)=50+50*(f1-f2)/(f1+f2);   %  The Frequency Slope
   %featvect(9)=100*halftime;    % time to half-peak-loudness
   %featvect(10)=cep_peakat;  % Location of the first non-trivial peak of the cepstrum
   featvect(6)=50+50*(sum1-sum2)/(sum1+sum2);  % Amplitude Slope;
   %featvect(12)=100*peakat;   % Location of peak amplitude
   %featvect(13)=10000*circ_mag; %Magnitude of harmonic-corrected frequency
   %featvect(14)=15*log10(spec_peak);    % Spectral peak                
   featvect(7)=100*spec_ent; %entropy of full spectrum


%   if any(imag(featvect))
%       tmp=find(imag(featvect));
%       disp('  Complex result in feature vector! ')
%       [tmp featvect(tmp)]
%   end
   
   return;


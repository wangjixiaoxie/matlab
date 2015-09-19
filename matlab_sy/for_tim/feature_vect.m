function featvect=feature_vect(syllable,Fs)

% featvect=FEATURE_VECT(syllable)
%
% Returns a row vector which is the feature vector of a particular syllable.
%

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

norm_flag=0;
norm_floor=90;  %number of decibels below max that is excluded from spect before norm

%values for calculating spectrogram
spect_win_dur = 10; %ms
spect_overlap = .50;  %percentage of overlap of specgram window
nfft=round(Fs*spect_win_dur/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points

% Calculate the spectral density of the syllable
[sdS,sdfS]=specgram(syllable,length(syllable),Fs,spect_win,noverlap);
sdFilt=sdS(find(sdfS<10000),:);    % Above 10000 Hz, filtered out anyway
sdf=sdfS(find(sdfS<10000));            %%%changed from 8k to match new filtering of sound files in autolabelkb
sd=abs(sdFilt);  % Throw away relative phase information
% Normalize
sd=sd/sum(sd);


%%%smooth waveform for mean amplitude calculation
smthsyl = conv(ones(10,1),abs(syllable));
offset = round((length(smthsyl)-length(syllable))/2);
smthsyl = smthsyl(offset:round(length(syllable))+offset);


% Calculate the spectral density of each half of the syllable
% [sds,sdfs]=specgram(syllable,1:ceil(length(syllable)),Fs,spect_win,noverlap);
% sd1=abs(sds(:,1));
% sd2=abs(sds(:,2));
% sum1=sum(sd1);
% sum2=sum(sd2);
% f1=sum(sd1.*sdfs)/sum1;
% f2=sum(sd2.*sdfs)/sum2;

%for pwr vs time calculation
lng1 = size(sdFilt,1)/2;
pf1=sum(abs(sdFilt(1:floor(lng1),:)));
pf2=sum(abs(sdFilt(floor(lng1)+1:2*lng1,:)));

%for pwr vs freq calculation
lng1 = size(sdFilt,2)/2;
pt1=sum(abs(sdFilt(:,1:floor(lng1))));
pt2=sum(abs(sdFilt(:,floor(lng1)+1:2*lng1)));

% Calculate the peak and spectral entropy of the full spectrogram
for i=1:size(sdS,1)*size(sdS,2)
    unwrap(i)=abs(sdS(i));
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
featvect(3)=10*log10(peak);                  % Third: Peak Loudness
featvect(4)=10*log2(Dur*1000);             % Fourth: Duration of Syllable (in ms);
featvect(5)=-100*sum(loud.*log2(loud))/log2(length(loud));    % Fifth: Entropy of the loudness vs. time (range: 0 to 100; 100 = flat amplitude; 0 = one very sharp peak)
featvect(6)=100*circ_avg;     %  A harmonic-corrected mean frequency (2^n+s^fv=fundamental frequency)
featvect(7)=100*circ_max;     %  A harmonic-corrected max frequency
featvect(8)=50+50*(sum(pf2)-sum(pf1))/(sum(pf1)+sum(pf2));   %  power vs frequency Slope
featvect(9)=100*halftime;    % time to half-peak-loudness
featvect(10)=cep_peakat;  % Location of the first non-trivial peak of the cepstrum
featvect(11)=50+50*(sum(pt2)-sum(pt1))/(sum(pt1)+sum(pt2));  % power vs time Slope;
featvect(12)=100*peakat;   % Location of peak amplitude
featvect(13)=10000*circ_mag; %Magnitude of harmonic-corrected frequency
featvect(14)=15*log10(spec_peak);    % Spectral peak
featvect(15)=100*spec_ent; %entropy of full spectrum
featvect(16)=10*log10(mean(smthsyl)); %mean amplitude


%   if any(imag(featvect))
%       tmp=find(imag(featvect));
%       disp('  Complex result in feature vector! ')
%       [tmp featvect(tmp)]
%   end

return;


function [sp t f]=mplotcbin(cbin,bnds,normfreqs)
%
% plots spectrogram of cbin file, which can be either .cbin or .wav
% also optionally returns spectrogram power (sp), time vector (t), and freq
% vector (f)
% 
% bnds = time limits, usually leave empty ([]) to plot whole file
% normfreqs isn't necessary
% example - mplotcbin('w39o13_021111_0657.-2158.cbin',[]);
%


[filepath,filename,fileext] = fileparts(cbin);

if(strcmpi(fileext,'.cbin'))
    [plainsong,fs] = ReadCbinFile(cbin);
elseif(strcmpi(fileext,'.wav'))
    [plainsong,fs] = wavread(cbin);
    plainsong = plainsong *10e3; % boost amplitude to cbin levels
end


% filter raw song
plainsong = bandpass(plainsong,fs,750,15000,'hanningfirff');

% provision for autoscaling if bounds are left unspecified. ~mnm 04.2009
if (isempty(bnds));
    bnds = [1/fs,(length(plainsong)*(1/fs))];
end

plainsong=ceil(plainsong(bnds(1)*fs:bnds(2)*fs));

%[sm,sp,t,f]=SmoothData(plainsong,fs);
%[sm,sp,t,f]=evsmooth(plainsong,fs,0,512,0.2,4); % SPTH(threshold), nFFT, degree FFT overlap, smooth window(ms)

[sp,t,f] = mFastSpect(plainsong,fs,0,256,0.9,2); % SPTH(threshold), nFFT, degree FFT overlap, smooth window(ms)

if(exist('normfreqs'));
    spt{1}=log(abs(sp));
    ps.freqlimbound=[6000 8000];
    ps.freqspacing=62.5;
    [colbnds]=calcbndsavn(spt,ps);
    [colbnds]=1.3*colbnds;
end

if(strcmpi(fileext,'.cbin'))
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
elseif(strcmpi(fileext,'.wav'))
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
end

if (exist('normfreqs'));
    set(gca,'CLim',colbnds);
end
hold on;

colormap('hot');
maxnorm=max(abs(plainsong));
% maxnorm=maxnorm/25;

%
% nrmsong=32768*(1/(1.01*maxnorm))*plainsong;
%
% %writeaiff_file
% aif_file='example.aiff'
% aiffwrite(aif_file,nrmsong,32000,16);
%
%
%  %  figure
% %  if (labels)
% %      load([fnvl{ii} '.not.mat'])
% %      zoomind=find(onsets/1000>bnds{ii}(1)&onsets/1000<bnds{ii}(2));
% %      for jj=1:length(zoomind)
% %          text(onsets(zoomind(jj))/1000-bnds{ii}(1), 200, labels(zoomind(jj)),'Color', 'w','Fontsize',16);
% % % ,'Color','w');
% %      end
%  end
%  box off;
%  axis on;
%
% outwav=wavwrite(nrmsong,32000, 16,'example.wav')

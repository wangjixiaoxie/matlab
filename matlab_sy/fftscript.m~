%This script will analyze and plot the raw spectra of wav files for entire
%song.




figure




wavfilestoread={'34.wav' '34m5.wav' '34m10.wav'}
colval={'r' 'b' 'y'}
init_time=2.15;
final_time=2.4;
%wavfilestoanalyze[2 6]


npts=4096;

figure

%this time for a restricted syllable
%clear dat;

%dat{1}=pshift7980{2};


%dat{1}=wavfilestoread{1}
%dat{2}=wavfilestoread{2};
%dat{3}=wavfilestoread{3};
%dat{4}=wavfilestoread{4};


loopval=1;
for i=[1 3 5]%]%th(wavfilestoread)
   
   dat{i}=sngshift{i};
    
  %[dat{i},fs,n]=wavread(wavfilestoread{i})
    fs=44100
    
    
    fftout{i}=fft(dat{i}(init_time*fs:final_time*fs))
    npts=length(fftout{i});
    pout{i}=fftout{i}.*conj(fftout{i})
    f = fs*(0:(npts/2))/npts;
    subplot(2,1,1)
    [sm,sp,t,fr]=evsmooth(dat{1}(init_time*fs:final_time*fs),fs,.05);
 imagesc(t,fr,log(abs(sp)));syn;ylim([0,1e4]);
  %  axis([3.2186 3.2798 0 7000])
    
    subplot(2,1,2)
    plot(f,(pout{i}(1:((npts/2)+1))), colval{loopval})
    loopval=loopval+1;
    hold on;
end
%legend(wavfilestoread{1},wavfilestoread{2},wavfilestoread{3},wavfilestoread{4},'TextColor','w')
%Now take the fft of each data file


%fftout{1}=fft(dat{1},npts)

%this gives power spectrum
%pout{1}=fftout{1}.*conj(fftout{1})/npts






%Graph the first 257 points (the other 255 points are redundant) on a meaningful frequency axis: 
%f = fs*(0:(npts/2))/npts;
%plot(f,log(pout{1}(1:((npts/2)+1))))
%title('Frequency content of y')



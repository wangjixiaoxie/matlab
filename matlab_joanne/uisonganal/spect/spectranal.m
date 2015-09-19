function [nspect, note_sn, best_freq, modulation, best_ratio, peak_vs_avg, f_centers]=spectranal(soundfile, norm_flag, filetype, d_flag)

% need to strip first bin of spectra: see psdanal for model that works
%Note: function deltaprods assumes that first bin of spectrum represents a frequency of f_step, NOT 0.



%normalization
%if 1st digit of norm_flag = 1 then, before any other normalization, divide each column of spect by avg spect.
%for second digit:
%normflag=x0; %no additional normalization
%norm_flag=x1;  %normalize to maximum value of each column
%norm_flag=x2;  %normalize to mean value of each column
%norm_flag=x3;  %normalize to average amplitude across all notes 
%norm_flag=x4;  %normalize to overall maximum


%minimum allowable signal to noise for notes that will be analyzed for noisiness (ratio of amplitudes)
min_sn=2;  %empirically this seems like a reasonable value

%minimum allowable value in spectrogram bin as percent of mean amplitude
% this can be set to zero or values above zero: if s/n is lousy, then can get infinite ratios
% (if many points of spectrum are zero after background subtraction)
% setting to small positive value will prevent this from happening, but also tend to make low amplitude notes "noisier"
% than high amplitude notes: not neccessarily innapropriate if done in conjunction with normalization to overall mean.
% need to think about/be aware of where this step is done

min_bin=.05; 

%values for filtering; filtering will be done after spectrogram: 
%uses a window that is symettric in frequency and has hanning shaped flanks
% is 0 at F_low hz and  F_high
% rises to a value of one over fband Hz  (see win_spect for more details);

F_low=500;
F_high=7000;
F_band=1500;

%values for calculating spectrogram
spect_win_dur=16;
spect_overlap = .80;  %percentage of overlap of specgram window
   
%paramaters for analysis: note that sinprods needs modification if fmin !=0
fmin_cut=0;  %don't change this without understanding what it does
fmax_cut=7000;

smooth_flag = 1;
gauss_wid = 50; %standard deviation of gauss filt (in Hz) used to smooth across freq

notefile=[soundfile,'.not.mat'];
filt_file=[soundfile,'.filt'];
spect_file=[soundfile,'.spect'];
   
   %get notefile
     if exist(notefile)   
       load(notefile);
     else
       disp(['cannot find ', notefile])
     end
 
%for now, don't use filt file, rather refilt with specified values    
  % get filt_file  
  % if exist(filt_file)     
   %  disp(['loading filtered song...',filt_file]);
    % [filtsong, Fs] = read_filt(filt_file);
     %   else
     default_Fs=Fs;
%     disp(['cannot find ', filt_file])
%     disp(['will attempt to filter songfile as default file type'])
     if nargin<=2; disp(['filetype is undefined']); end;
     [rawsong,Fs]=soundin([],soundfile, filetype);
%     disp('filtering song...');
     %unless Fs has been read from the file use the default value which was
     %set by the user or previously read in from a note file
     if Fs == -1; Fs = default_Fs; end
     %filter soundfile
 %    filtsong=bandpass(rawsong,Fs,F_low,F_high);
 %  end
 
 filtsong=rawsong;
  
   
   %calculate spectrogram
    disp('calculating spectrogram...');
    %first calculate nfft and noverlap
     nfft=round(Fs*spect_win_dur/1000);
     spect_win = hanning(nfft);
     noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points       
    %now calculate spectrogram
     [spect, freq, time] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
     spect=abs(spect);
     %idx_spect=scale_spect(spect);  %calculate index array for spectrogram
     f_min = freq(1);
     f_max = freq(length(freq));
     freq_spect = [f_min, f_max];
     t_min = time(1)*1000; %convert to ms
     t_max = time(length(time))*1000; %convert to ms
    %adjust time axis for spectrogram offset (1/2 window duration in ms)
     t_min = t_min + 0.5*nfft*1000/Fs;  
     t_max = t_max + 0.5*nfft*1000/Fs;  
     time_spect = [t_min, t_max];  
   nfbins=size(spect,1);
   f_step=(f_max-f_min)/(nfbins-1);
     
   [spect]=win_spect(spect,f_step,F_low,F_high,F_band, 1);          
  
   %separate spectrogram file into background and notes
   % within program split_spect have option of constraining spectra to completely within note or only requiring 
   % center of window to fall withinnote.
   [note_spect, back]=split_spect(spect,t_min,t_max, noverlap, spect_win, onsets, offsets);
   
   
   %calculate mean background from quietest 50% of background (because of flapping,etc.)
   med_back=median(mean(back));
   i=find(mean(back)<med_back);
   back=back(:,i);
   back=mean(back');
   back=back';

   nfbins=size(spect,1);
   f_step=(f_max-f_min)/(nfbins-1);

  %smooth across frequencies
  if smooth_flag == 1
      %calculate gaussian filter out to +/- 3std dev at appropriate intervals
      nvals=ceil(3*gauss_wid/f_step);
      interp_vals = [-1*f_step*nvals:f_step:f_step*nvals];
      gauss_filt = normpdf(interp_vals,0,gauss_wid);
      gauss_filt = gauss_filt/sum(gauss_filt);
      gauss_filt = gauss_filt';
      edge_trash = (length(gauss_filt)-1)/2;
      sback=conv(gauss_filt,back);
      back=sback(edge_trash+1:length(sback)- edge_trash);
 %*****************************************************
 %  next line inserted to return whole spectrogram instead of just note_spect
 %******************************************************     
      note_spect = spect;
      note_spect = conv2(note_spect,gauss_filt,'same'); 
   end
   
     
 
   %cut to desired range of frequencies
   f_centers=[0:nfbins-1]*f_step+f_min;
   keep_f=find(f_centers>=fmin_cut & f_centers <= fmax_cut);
   back=back(keep_f,:);
   note_spect=note_spect(keep_f,:);
   f_centers=f_centers(keep_f);

   %calculate average signal/noise ratio for each point in note_spect.
   mean_noise=mean(back);
   mean_signal = mean(note_spect);
   note_sn=mean_signal/mean_noise;
  
   %eliminate notes (spectra) that do not meet appropriate signal to noise criterion
   [keep_idx]=find(note_sn>=min_sn);
   note_spect=note_spect(:,keep_idx);
   mean_signal=mean_signal(keep_idx);
   note_sn=note_sn(keep_idx);
    
   %subtract background from notes
   sback=back*ones(1,size(note_spect,2));
   snote_spect=note_spect-sback;
   
   %set noise floor to zero
   snote_spect=max(0,snote_spect); 
   
if norm_flag >= 10;
    %divide the note spectrum by the mean note spectrum
    mspect=mean(snote_spect');
    %don't divide by zero
    mspect=max(.01*mean(mspect),mspect);
    norm_mat=mspect'*ones(1,size(snote_spect,2));
    snote_spect=snote_spect./norm_mat;
end

%refilter here because of potentially large values at edges due to ratio
[spect]=win_spect(spect,f_step,F_low,F_high,F_band, 1);        

if norm_flag == 10;
  nspect=snote_spect;
end
           



 

 if norm_flag == 1 | norm_flag == 11;
   %normalize each column of note_spect to the maximum intensity for that column
   max_vals=max(snote_spect);
   norm_fact=1./max_vals;
   norm_diag=diag(norm_fact);
   nspect=snote_spect*norm_diag;
 end
   
 if norm_flag == 2 | norm_flag == 12;  
   %normalize each column of note_spect to the avg amplitude for that column
   mean_vals=mean(snote_spect);
   norm_fact=1./mean_vals;
   norm_diag=diag(norm_fact);
   nspect=snote_spect*norm_diag;
 end
   
 if norm_flag == 3 | norm_flag == 13; 
   %normalize by average amplitude across all notes
   meanamp=mean(mean(snote_spect));
   nspect=snote_spect/meanamp;
 end
   
  if norm_flag == 4 | norm_flag == 14; 
   %normalize by greatest amplitude across all notes
   maxamp=max(mean(snote_spect));
   nspect=snote_spect/maxamp;
 end
     

   %set floor to minimum allowable value
   overall_mean=mean(mean(nspect));
   bin_floor=min_bin*overall_mean;
   nspect=max(bin_floor,nspect);


   %measure modulation for each column
   modulation=sum(abs(diff(nspect)));
   modulation=modulation/(max(f_centers)-f_centers(2));
   
   %measure best sinprod for each column

   %[cos_matrix, cos_prods, best_ratio, best_per, peak_vs_avg]=deltaprods(nspect,f_step);
    
   %best_freq=best_per.*f_step;
 

  
if d_flag ==1  
   %disply some stuff
  
   figure
   subplot(3,1,1)
   imagesc(flipud(snote_spect));
   hold on;
   bins=size(nspect,1);
   plotbf=bins-best_per;
   plot(plotbf,'b');
   hold off;
   subplot(3,1,2)
   plot(peak_vs_avg);
   ylabel('peak_vs_avg');
   axis([0, length(peak_vs_avg), 0, max(peak_vs_avg')]);
   subplot(3,1,3);
   plot(log(best_ratio));
   ylabel('log(best ratio)');
   axis([0, length(best_ratio), 0, log(max(best_ratio'))]);
end

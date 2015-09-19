function [best_freq, best_aratio,best_pratio, ampspect, f_centers]=psdanal(song,Fs);
%function [best_freq, best_aratio,best_pratio, modulation, ampspect, f_centers]=psdanal(song,Fs);

%psdanal takes segment of song and sample rate and displays psd as
% well as other stuff such as estimate of ff.  Much of this
% will be pirated from spectranal
% for now, no background subtraction and fixed normalization, but could add





global percent_from_start ms_from_start psd_label soundfile sampled_dur

%parameters to be used for filtering and display:

%filter spectrum because we have edge splatter in chopped input signal
%values for filtering; filtering will be done to spectrum: 
%uses a window that is symettric in frequency and has hanning shaped flanks
% is 0 at F_low hz and  F_high
% rises to a value of one over fband Hz  (see win_spect for more details);

F_low=300;
F_high=8100;
F_band=200;


%values for calculating spectrum
spect_win_dur=16;
spect_overlap = .50;  %percentage of overlap of specgram window
   
%paramaters for analysis: note that deltaprods needs modification if fmin !=0
fmin_cut=0;  %don't change this without understanding what it does
fmax_cut=7000;

%if log_flag==1 then use log of power = amp
log_flag=1;


smooth_flag = 1;
gauss_wid = 50; %standard deviation of gauss filt (in Hz) used to smooth across freq



 
   %calculate spectrum
    disp('calculating spectrum...');
    %first calculate nfft and noverlap
     nfft=round(Fs*spect_win_dur/1000);
     spect_win = hanning(nfft);
     noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points       
    %now calculate spectrum
     [spect,pcxx, freq] = psd(song,nfft,Fs,spect_win,noverlap,[]);
     f_min = freq(1);
     f_max = freq(length(freq));
     freq_spect = [f_min, f_max];
     nfbins=size(spect,1);
     f_step=(f_max-f_min)/(nfbins-1);
     
     
     %
     %filter spectra
     [spect]=win_spect(spect,f_step,F_low,F_high,F_band, 1);
     
     
     
     %normalize power by dividing psd by mean psd.
     avg_power=mean(spect);
     spect=spect/avg_power;
     
         
    %smooth across frequencies
      if smooth_flag == 1
        %calculate gaussian filter out to +/- 3std dev at appropriate intervals
        nvals=ceil(3*gauss_wid/f_step);
        interp_vals = [-1*f_step*nvals:f_step:f_step*nvals];
        gauss_filt = normpdf(interp_vals,0,gauss_wid);
        gauss_filt = gauss_filt/sum(gauss_filt);
        gauss_filt = gauss_filt';
        edge_trash = (length(gauss_filt)-1)/2;
        %sback=conv(gauss_filt,back);
        %back=sback(edge_trash+1:length(sback)- edge_trash);
        spect = conv2(spect,gauss_filt,'same'); 
      end
      
    %cut to desired range of frequencies for analysis
        f_centers=[0:nfbins-1]*f_step+f_min;
        keep_f=find(f_centers>=fmin_cut & f_centers <= fmax_cut);
        %back=back(keep_f,:);
       spect=spect(keep_f,:);
        f_centers=f_centers(keep_f);
    
    if log_flag==1; 
     %convert psd to amplitude 
     ampspect=spect;
     %set min val prior to taking log
     zero_ind=find(ampspect==0);
     non_zero_ind=find(ampspect>0);
     min_val=min(ampspect(non_zero_ind));
     ampspect(zero_ind)=min_val*ones(size(zero_ind));
     ampspect=log(ampspect);
    end
     
     %measure modulation
     %   modulation=sum(abs(diff(ampspect)));
      %  modulation=modulation/(max(f_centers)-f_centers(2));
   
     %measure best deltaprod
        %need to strip zero f bin of spectrum for deltaprods to work correctly
        nfbins=size(ampspect,1);
	short_ampspect=ampspect(2:nfbins,:);
        [best_freq, best_idx, ratio_prods, ratio_freq_scale]=deltaprods(short_ampspect,f_step);
  diff_prods=ratio_prods(:,1);
  percent_harm=ratio_prods(:,2);
  best_aratio=diff_prods(best_idx);
  best_pratio=percent_harm(best_idx);
  %plot stuff
  h_psdwin=figure;
  h_ratio_plot=subplot(2,1,1);
  h_r_data=plot(ratio_freq_scale,diff_prods,'g');
  set(gca,'ylim',[0,1.2*max(diff_prods)]);
  hold on
  h_p_data=plot(ratio_freq_scale,percent_harm,'m');
  h_f_line=plot([best_freq, best_freq],[0, best_aratio],'r');
  title(['best aratio = ',num2str(best_aratio)]);
  h_spect_plot=subplot(2,1,2);
  plot(f_centers,ampspect,'g');
  hold on;
  %pspect=spect*max(ampspect)/max(spect);
  %plot(f_centers,pspect,'r');
  n=0;
  for i = best_freq:best_freq:max(f_centers)
    n=n+1;
    h_f_lines(n)=plot([i,i],[min(ampspect),max(ampspect)],'r');
  end
  title(['best frq = ',num2str(best_freq)]);
  %save axes handles in window userdata
  set(h_psdwin,'userdata',[h_ratio_plot,h_spect_plot]);
  %save fline handles with apprpriate axes userdata
  set(h_ratio_plot,'userdata',[h_f_line, h_r_data,h_p_data]);
  set(h_spect_plot,'userdata',[h_f_lines]);
  
  %set ratio_prods axis button down function to allow user changing of F
  set(h_ratio_plot,'buttondownfcn',[...
               'p=get(gca,''currentpoint'');'...
	       'best_freq=p(1,1);',...
	       'psdupdate(gcf,best_freq);']);
	      
	      %'handles=get(gcf,''userdata'');'...
	      %'subplot(handles(1));'...
	      %'h_f_line=get(gca,''userdata'');'...
	      %'set(h_f_line,''xdata'',[best_freq,best_freq]);'...
	      %'])

  %display results to screen (or file)
  disp('soundfile   note  bfreq  a_ratio p_ratio  %in ms_in ms_dur')
  fprintf(1,'\n');
  fprintf(1,'%s\t',[soundfile]);
  fprintf(1,'%s\t',[psd_label]);
  fprintf(1,'%5.2f\t',[best_freq, best_aratio, best_pratio]);
  fprintf(1,'%4.1f\t',[percent_from_start, ms_from_start, sampled_dur]);  
  fprintf(1,'\n');










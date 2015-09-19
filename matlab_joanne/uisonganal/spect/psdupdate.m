function psdupdate(h_figure,best_freq)

%this function is called by clicking on the best_ratios window of psdanal plot
% it updates the frequency display and  prints values to screen or file...
% It is assumed that figure userdata has handles to subplots
%and that subplot userdata has handles to frequency lines.

global percent_from_start ms_from_start psd_label soundfile sampled_dur

%get handles
  handles=get(h_figure,'userdata');
  h_ratio_plot=handles(1);
  h_spect_plot=handles(2);
  
  ratio_handles=get(h_ratio_plot,'userdata');
  h_f_line=ratio_handles(1);
  h_r_data=ratio_handles(2);
  h_p_data=ratio_handles(3);
  h_f_lines=get(h_spect_plot,'userdata');
  
  
%move freq_line on ratio plot
 set(h_f_line,'xdata',[best_freq,best_freq]);
 
%get new ratio vals
 subplot(h_ratio_plot);
 diff_prods=get(h_r_data,'ydata');
 percent_harm=get(h_p_data,'ydata');
 ratio_freq_vals=get(h_r_data,'xdata'); 
 %find closest freq bin
 [y, bf_ind]=min(abs(ratio_freq_vals-best_freq));
 best_aratio=diff_prods(bf_ind);
 best_pratio=percent_harm(bf_ind);
 title(['best aratio = ',num2str(best_aratio)]);

%move frequency lines on spect_plot
  subplot(h_spect_plot);
  xl=get(gca,'xlim');
  max_freq=xl(2);
  yl=get(gca,'ylim');
  n=0;
  for i = best_freq:best_freq:max_freq
    n=n+1;
    h_f_lines1(n)=plot([i,i],[yl(1),yl(2)],'r');
  end
  delete(h_f_lines);
  set(h_spect_plot,'userdata',h_f_lines1);
  title(['best frq = ',num2str(best_freq)]);


  
  
%display data to screen

%display results to screen (or file)
  

  fprintf(1,'\n');
  fprintf(1,'%s\t',[soundfile]);
  fprintf(1,'%s\t',[psd_label]);
  fprintf(1,'%5.2f\t',[best_freq, best_aratio, best_pratio]);
  fprintf(1,'%4.1f\t',[percent_from_start, ms_from_start, sampled_dur]);  
  fprintf(1,'\n');
  

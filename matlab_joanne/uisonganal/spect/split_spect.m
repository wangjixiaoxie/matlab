function [note_spect, back_spect]=split_spect(idx_spect,t_min,t_max, noverlap, spect_win, onsets, offsets)


%should bins be entirely within notes?
strict_flag = 1;


[nfreqs,nbins] = size(idx_spect);
step_size=(t_max-t_min)/(nbins-1);
bin_centers=[0:nbins-1]*step_size+t_min;
percent_overlap=noverlap/length(spect_win);
bin_dur=step_size/(1-percent_overlap);
bin_starts=bin_centers-.5*bin_dur;
bin_ends=bin_centers+.5*bin_dur;

if strict_flag==1
 %get subset of spectrogram that contains bins that are entirely within notes

 for i=1:length(onsets)
  after_on_ind=find(bin_starts>onsets(i));
  before_off_ind=find(bin_ends<offsets(i));
  startnote=after_on_ind(1);
  endnote=max(before_off_ind);
  note_ind=[startnote:endnote];
  curr_note=idx_spect(:,note_ind);
  note_spect=[note_spect,curr_note];
 end

 onsets=[onsets; t_max+bin_dur];
 %get background
  %first interval
  before_on=find(bin_ends<onsets(1));
  back_spect=idx_spect(:,1:max(before_on));
 for i=1:length(offsets)
  after_off_ind=find(bin_starts>offsets(i));
  before_on_ind=find(bin_ends<onsets(i+1));
  startint=after_off_ind(1);
  endint=max(before_on_ind);
  int_ind=[startint:endint];
  curr_int=idx_spect(:,int_ind);
  back_spect=[back_spect,curr_int];
 end
else
 %get subset of spectrogram that has bins with centers contained in notes
 
 for i=1:length(onsets)
  after_on_ind=find(bin_centers>onsets(i));
  before_off_ind=find(bin_centers<offsets(i));
  startnote=after_on_ind(1);
  endnote=max(before_off_ind);
  note_ind=[startnote:endnote];
  curr_note=idx_spect(:,note_ind);
  note_spect=[note_spect,curr_note];
 end

 %still use strict criterea for intervals
 onsets=[onsets; t_max+bin_dur];
 %get background
  %first interval
  before_on=find(bin_ends<onsets(1));
  back_spect=idx_spect(:,1:max(before_on));
 for i=1:length(offsets)
  after_off_ind=find(bin_starts>offsets(i));
  before_on_ind=find(bin_ends<onsets(i+1));
  startint=after_off_ind(1);
  endint=max(before_on_ind);
  int_ind=[startint:endint];
  curr_int=idx_spect(:,int_ind);
  back_spect=[back_spect,curr_int];
 end

end




  

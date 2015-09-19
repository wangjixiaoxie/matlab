function disp_idx_spect(idx_spect, time_spect, freq_spect, spect_floor, spect_ceil, spect_range, varargin)

%display indexed representation of spectrogram (as creatd by scale_spect and saved in spect files)
% the times and frequencies associated with displayed image are set by 
% time_spect = [t_min t_max];   freq_spect = [f_min f_max]
% the spectrogram is assumed to be in indexed representation such that values within n dB of max have an index of n
% the spectrogram floor is asumed to be -100 dB 
%floor_db set the minimum value displayed in db relative to the maximum value
%ceil_db sets the value at which display saturates in in db relative to the maximum value 


%build appropriate colormap
if nargin == 6
  cm = make_map(spect_floor, spect_ceil, spect_range);
elseif nargin == 7
  spect_cmap_name = varargin{1};
  cm = make_map(spect_floor, spect_ceil, spect_range, spect_cmap_name);
elseif nargin == 8
  spect_cmap_name = varargin{1};
  spect_gamma_type = varargin{2};
  cm = make_map(spect_floor, spect_ceil, spect_range, ...
      spect_cmap_name, spect_gamma_type);
end

%display that sucker
%because of matlab bug, cannot zoom correctly on image that has axis = xy
% get message about warped color image, and boundaries of image columns jump around
%as work around, invert image data and display with images default reversd y axis
%this has the consequence of screwing up y axis lables, which plot from smallest to largest along the y axis direction
%funky workaround is to display frequency with negative numbers
%idx_spect = flipud(idx_spect);
%freq_spect = fliplr(freq_spect*(-1));

cla
hold on
image(time_spect, freq_spect, idx_spect), colormap(cm);
hold off


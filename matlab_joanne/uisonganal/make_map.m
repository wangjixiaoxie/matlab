function spect_map = make_map(spect_floor, spect_ceil, spect_range, varargin)

%build colormap for display of idx_spect
%spectrogram indices indicate number of dB from max value
% see scale_spect and disp_idx_spect for more info

% Note: 'classic' gamma_type only valid for grayscale colormap
if nargin >= 4
  cmap_name = varargin{1};
  spect_gamma_type = 'brightness';
elseif nargin == 3
  % Use default colormap
  spect_gamma_type = 'classic';
end
if nargin == 5
  if strncmpi(cmap_name,'gray',4)    
    spect_gamma_type = varargin{2};
  end
end
  
spect_floor = round(abs(spect_floor));
spect_ceil = round(abs(spect_ceil));

%default in case of screwy values is spect_ceil = 1; use full color range
if spect_ceil < 1 | spect_ceil > 100
   spect_ceil = 1;
end

%spect_floor must be greater than or equal to spect_ceil
if spect_floor < spect_ceil
   spect_floor = spect_ceil;
end 

%spect range must be between .01 and 1
%spect_range = min(max(spect_range, .01), 1);

%build appropriate colormap

%indexes 1:spect_ceil will result in saturation
saturated = zeros(spect_ceil-1,3); 

%indexes spect_floor:100 will be blanked
blanked = ones(99-spect_floor,3);

%continual gradation of gray scale for indeces in between
map_size = 100 - size(saturated,1) - size(blanked,1);

if strncmpi(spect_gamma_type,'classic',7)
  %max_val = spect_range;
  %gray1 = [max_val/map_size:max_val/map_size:max_val]';
  %min_val =  1-max_val;
  gray1 = flipud(1-([1/map_size:1/map_size:1].^spect_range)');
  cmap = [gray1,gray1,gray1];  
else
  % Map beta between -1 and 1.
  % Note: spect_range is between 0 and 4, with default at 2.
  beta = (spect_range - 2)/2;
  eval(['cmap = ', cmap_name, '(', num2str(map_size), ');']);
  cmap = brighten(cmap,beta);
end

spect_map = ([saturated; cmap; blanked]);

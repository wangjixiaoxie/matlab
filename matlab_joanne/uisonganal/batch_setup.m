function [interact] = batch_setup

%gets parameters for starting batch processing

%control flow variables that are set by batch setup
global save_notes save_filt save_spect batch_disp batch_print do_filt
                  
%set options
disp('')
disp('Set options for non-interactive batch processing')
disp('')

%save notefile?
if save_notes == 1
  temp = 'y';
else
  temp = 'n';
end
disp([' Save notefile? [', temp, ']'])
temp = input('','s');
if strcmp(temp,'y') | strcmp(temp,'Y') | strcmp(temp,'yes')
  save_notes = 1;
elseif strcmp(temp,'n') | strcmp(temp,'N') | strcmp(temp,'no')
  save_notes = 0;                  
end
         

%do filtering?
if do_filt == 1
  temp = 'y';
else
  temp = 'n';
end
disp([' Filter songs? [', temp, ']'])
temp = input('','s');
if strncmpi(temp,'y',1)
  do_filt = 1;
elseif strncmpi(temp,'n',1)
  do_filt = 0;                  
end
         
%save filtered song?
if do_filt
  if save_filt == 1
    temp = 'y';
  else
    temp = 'n';
  end
  disp([' Save filtered song? [', temp, ']'])
  temp = input('','s');
  if strcmp(temp,'y') | strcmp(temp,'Y') | strcmp(temp,'yes')
    save_filt = 1;
  elseif strcmp(temp,'n') | strcmp(temp,'N') | strcmp(temp,'no')
    save_filt = 0;                  
  end
end

%save spectrogram?
if save_spect == 1
  temp = 'y';
else
  temp = 'n';
end
disp([' Save spectrogram? [', temp, ']'])
temp = input('','s');
if strcmp(temp,'y') | strcmp(temp,'Y') | strcmp(temp,'yes')
  save_spect = 1;
elseif strcmp(temp,'n') | strcmp(temp,'N') | strcmp(temp,'no')
  save_spect = 0;                  
end
         
%display while running?
if batch_disp == 1
  temp = 'y';
else
  temp = 'n';
end
disp([' Display data while running? [', temp, ']'])
temp = input('','s');
if strcmp(temp,'y') | strcmp(temp,'Y') | strcmp(temp,'yes')
  batch_disp = 1;
elseif strcmp(temp,'n') | strcmp(temp,'N') | strcmp(temp,'no')
  batch_disp = 0;                  
end

%print spectrograms?
if batch_print == 1
  temp = 'y';
else
  temp = 'n';
end
disp([' Print spectrogram/amplitude plots? [', temp, ']'])
temp = input('','s');
if strcmp(temp,'y') | strcmp(temp,'Y') | strcmp(temp,'yes')
  batch_print = 1;
elseif strcmp(temp,'n') | strcmp(temp,'N') | strcmp(temp,'no')
  batch_print = 0;                  
end         

%start batch processing
disp(' Start batch processing?  [y]')
temp = input('','s');
if strcmp(temp,'y') | strcmp(temp,'Y') | strcmp(temp,'yes') | isempty(temp)
  interact = 0;
elseif strcmp(temp,'n') | strcmp(temp,'N') | strcmp(temp,'no') | strcmp(temp,'q') | strcmp(temp,'Q')
  interact = 1;         
end
     
     

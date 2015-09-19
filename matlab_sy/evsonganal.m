function varargout = evsonganal(varargin)

% set up sub paths
mypath = fileparts(which('uisonganal')) % This might not be foolproof!
subpath = fullfile(mypath,'spect','')
addpath(subpath)

%some global variables   

%values for display
global d_song time_d_song     %values used by disp_song
global spect time_spect freq_spect
global threshold smooth Fs
global h_main h_main_spect h_main_amp h_main_labels
global win_percent
global pr_left_margin pr_bottom_margin
global main_win_code
global Fs
global spect_floor spect_ceil spect_range spect_cmap_name spect_gamma_type
global soundfile
global do_zc d_zc
global zc zcwin

%control flow
global get_next interact batch_disp save_notes save_spect ...
    save_filt do_filt batch_print

%values used in segmentation
global threshold min_int min_dur sm_win
global onsets offsets labels
global smooth filtsong

%file data and pathnames
global soundfile path_songfile path_notefile path_filtfile path_spectfile

%indexes into main figure userdata
%these just facilitate mnemonics of storage and retrieval of values in userdata
%values are set below and used by programs that retrieve from userdata
global imain_win_code idata_status ih_main_amp ih_main_spect ih_main_labels ih_fname 
global ih_path_songfile ih_path_notefile ih_path_filtfile ih_path_spectfile ih_amp_plot ih_zc_plot 

%some default values

%identification code for figure containing spectrogram and amplitude plots
% stored as first value of figure userdata
main_win_code = -101;

%indexes for userdata storage and retrieval
%code identifying window type (number)
imain_win_code = 1;
%code indicating status of data for current figure
%see make_current for meaning of status codes
idata_status = 2;
%handles to amp and spect axes
ih_main_amp = 3;
ih_main_spect = 4;
%handles to hidden text containing filename and path strings
ih_fname = 5;
ih_path_songfile = 6;
ih_path_notefile = 7;
ih_path_filtfile = 8;
ih_path_spectfile = 9;
%handles to amplitude plot
ih_amp_plot = 10;
% for the label bar
ih_main_labels = 11;
%handles to the zero crossings plot
ih_zc_plot = 12;

% Do we use some new gui features?
use_new_gui_features = 1

%main window, starts invisible
% In an axis tagged "SAAxis", set the pointer to a cross
h_main=figure('Tag','SAMain','WindowButtonMotionFcn','centercb(''center_pointer'')','WindowButtonDownFcn','centercb(''center'')', 'Visible', 'off');
%h_main = [];
varargout{1} = h_main;
SA.handles.SAMain = h_main;
% Flag which determines whether to get the next song.
get_next = 0;

%*******************************************
% Default values for user parameters
%*******************************************
% To change parameters, save new preferences rather than editing this file!

%interactive or batch mode?
SA.param_defaults.interact = 1;

%display data during batch processing?
SA.param_defaults.batch_disp = 1;

%print in batch mode?
SA.param_defaults.batch_print = 0;

%save notefiles, filtered song, spectrogram?
SA.param_defaults.save_notes = 1;
SA.param_defaults.save_filt = 1;
SA.param_defaults.save_spect = 1;

%do we need to filter the data?
SA.param_defaults.do_filt = 1;

%for bandpass
SA.param_defaults.F_low=300;
SA.param_defaults.F_high=8000;
SA.param_defaults.filter_type = 'hanning';

%window length for smoothing prior to segmentation (in ms)
SA.param_defaults.sm_win = 2.0;

%threshold for segmentation, a value of 0 causes threshold to be calculated
SA.param_defaults.threshold = 0;     %is there a current threshold? 

% Do we display zero crossing frequency?
SA.param_defaults.do_zc = 0;

% Window used for zero crossing estimation (in ms)
SA.param_defaults.zcwin = 10.0;

%d_Fs is resample rate for more rapid display of song
SA.param_defaults.d_Fs=1000;

%criteria for intervals and notes (minimum durations in ms)
SA.param_defaults.min_int=5;
SA.param_defaults.min_dur=20;
%SA.param_defaults.max_dur=120;

%values for calculating spectrogram
SA.param_defaults.spect_win_dur=16;
SA.param_defaults.spect_overlap = .80;  %percentage of overlap of specgram window

%initial floor, ceiling and contrast for display of spectrogram
SA.param_defaults.spect_floor = -55;  %floor cutoff in dB re max
SA.param_defaults.spect_ceil = 0;   %ceil saturation in db re max
% Note: spect_range is a between 0 and 4, with default at 2.
% For brightness mapping the beta in brighten(cmap,beta) is
% (spect_range - 2)/2 (i.e., from -1:1);
SA.param_defaults.spect_range = 2.0; 
% Classic mode uses a funky sort of gamma correction. Other option
% is 'brightness'.
SA.param_defaults.spect_gamma_type = 'classic';
% Non grayscale colormaps can only use the 'brightness' gamma_type.
SA.param_defaults.spect_cmap_name = 'gray';

% Default range on the time axis in msec. (-1) means full range;
SA.param_defaults.xrange_start = -1;

%win_percent sets the amount by which window is shifted left or right 
win_percent = 0.80;

%default values for printing
SA.param_defaults.pr_orient = 'landscape';  %paper orientation
SA.param_defaults.pr_papertype = 'usletter'; %papertype 
SA.param_defaults.pr_y_paperinches = 6.5; %number of inches for height of printed figure
SA.param_defaults.pr_x_per_paperinch = 600; %number of units(ms) per inch of printed figure
SA.param_defaults.pr_left_margin = 1.0;
SA.param_defaults.pr_bottom_margin = 1.0;

% End of user parameter defaults

SA.samprate = 32000;
SA.params_loaded = 0;
SA.file_loaded = 0;
setappdata(h_main, 'SA_Data', SA);
       
if use_new_gui_features
  if nargin == 1
    config_file = varargin{1};
  end

  % New gui setup (make sure these are modal dialogs for now)

  % Load file dialog
  filecb('loadfile', h_main);
  SA = getappdata(h_main, 'SA_Data');

  % Do we continue?
  if ~SA.file_loaded
    error('No file loaded! Goodbye.')
  end
  % Check directories for sanity. 
  % Set a sane default if found to be invalid.
  [SA.loaddir, fname,fext,fver] = fileparts(SA.loadfile);
  if ~isfield(SA,'resultsdir') | isempty(SA.resultsdir) | ...
	~exist(SA.resultsdir,'dir') 
    SA.resultsdir = SA.loaddir
  end
  if strcmp(SA.resultsdir,'.')
    SA.resultsdir = pwd;
  end
  if ~isfield(SA,'tempdir') | isempty(SA.tempdir) | ...
	~exist(SA.tempdir,'dir') 
    SA.tempdir = SA.resultsdir
  end
  if strcmp(SA.tempdir,'.')
    SA.tempdir = pwd;
  end
  if ~isfield(SA,'rawdatadir') | isempty(SA.rawdatadir) | ...
	~exist(SA.rawdatadir,'dir') 
    SA.rawdatadir = SA.loaddir
  end
  if strcmp(SA.rawdatadir,'.')
    SA.rawdatadir = pwd;
  end  
  setappdata(h_main, 'SA_Data', SA);
    
  % Bring up parameter setting dialog
  if exist('config_file','var') & exist(config_file,'file')
    paramcb('construct',h_main,config_file);
  else
    paramcb('construct',h_main);
  end

  SA = getappdata(h_main, 'SA_Data');
    
  params = fieldnames(SA.param_defaults)
  if SA.params_loaded
    for ip=1:length(params)
      % For testing and compatibility with old way
      eval([params{ip}, ' = SA.params.', params{ip}])
    end
  else
    % Load in defaults
    for ip=1:length(params)
      % For testing and compatibility with old way
      eval([params{ip}, ' = SA.param_defaults.', params{ip}])
    end
  end
  
else 
  %get name of metafile containing songfile names

  SA.single = 0;
  metafile = 0;
  while metafile == 0 | isempty(metafile)
    disp('select batchfile');
    [metafile, SA.loaddir]=uigetfile('*','select batchfile')
    SA.loadfile = fullfile(SA.loaddir,metafile)
  end
  
  % set paths for storage and retrieval, for now default to pathname,
  % that is, the same place where the batch file is, but
  % be more flexible here and let things be in different places.
  
  disp(['You must enter the directory for uisonganal results and temporary files below;'])
  disp(['to use the current directory enter ".";'])
  SA.resultsdir = input(['Results directory [', SA.loaddir, ']: '],'s')
  if isempty(SA.resultsdir)  
    SA.resultsdir = SA.loaddir;
  elseif strcmp(SA.resultsdir,'.')
    SA.resultsdir = pwd;
  end
  % Check for valid path
  while ~exist(SA.resultsdir,'dir')
    disp('Directory does not exist. Try again.')
    SA.resultsdir = input(['Results directory [', SA.loaddir, ']: '],'s')
    if isempty(SA.resultsdir)  
      SA.resultsdir = SA.loaddir;
    elseif strcmp(SA.resultsdir,'.')
      SA.resultsdir = pwd;
    end
  end  
  SA.tempdir = SA.resultsdir
  
  disp(['You must enter the base directory for song files below;'])
  disp(['If the batch file has absolute paths in it enter "/";'])
  disp(['to use the current directory enter ".";'])
  SA.rawdatadir = input(['Song file base directory [', SA.loaddir, ']: '],'s')
  if isempty(SA.rawdatadir)  
    SA.rawdatadir = SA.loaddir
  elseif strcmp(SA.rawdatadir,'/')
    SA.rawdatadir = ''
  elseif strcmp(SA.rawdatadir,'.')
    SA.rawdatadir = pwd  
  else
    % Check for valid path
    while ~exist(SA.rawdatadir,'dir')
      disp('Directory does not exist. Try again.')
      SA.rawdatadir = input(['Song file base directory [', SA.loaddir, ']: '],'s')
      if isempty(SA.rawdatadir)  
	SA.rawdatadir = SA.loaddir
      elseif strcmp(SA.rawdatadir,'/')
	SA.rawdatadir = ''
	break
      elseif strcmp(SA.rawdatadir,'.')
	SA.rawdatadir = pwd  
      end
    end
  end
  
  %get filetype, later versions of read program can be smarter about this
  
  disp('What is type of sound file? [w]')
  disp(' b = binary from mac')
  disp(' w = wavefile (i.e. cbs)')
  disp(' d = dcpfile')
  disp(' filt = .filt file from uisonganal')
  disp(' foo = foosong/gogo (sets Fs = 32000)')
  disp(' o = observer file (last/song channel)')
  disp(' o1r = observer file (second to last)')
  
  SA.filetype = 'null';
  while strcmp(SA.filetype,'null')
    temp=input(' ','s');
    if strncmpi(temp,'b',1);
      SA.filetype = 'b';
    elseif strncmpi(temp,'w',1) | isempty(temp)
      SA.filetype = 'w';
    elseif strncmpi(temp,'d',1)
      SA.filetype = 'd';  
    elseif strncmpi(temp,'filt',4)
      SA.filetype = 'filt';  
    elseif strncmpi(temp,'foo',3)
      SA.filetype = 'foo';  
    elseif strcmpi(temp,'o')
      SA.filetype = 'obs0r';
    elseif strcmpi(temp,'o1r')
      SA.filetype = 'obs1r';
    else
      disp('Unacceptable! Pick again')
      SA.filetype = 'null';
    end
  end
  
  %get default sample rate
  % sampling rate of WAV files gets loaded in from the file automatically 
  SA.samprate = 32000;
  if ~(strcmp(SA.filetype,'w') | strcmp(SA.filetype,'filt'))
    SA.samprate = input('What is sample rate (samples/sec)? [32000] ');
    if isempty(SA.samprate); SA.samprate = 32000; end
  end

  % Set all parameters to defaults
  params = fieldnames(SA.param_defaults)
  for ip=1:length(params)
    eval(['SA.params.', params{ip}, ' = SA.param_defaults.', params{ip}])
    % For testing and compatibility with old way
    eval([params{ip}, ' = SA.param_defaults.', params{ip}])
  end

  % Save the settings in the UI's appdata
  setappdata(h_main, 'SA_Data', SA);

end % if use_new_gui_features

% Always display in interactive mode
if interact
  batch_disp = 1;
  SA.batch_disp = batch_disp;
end

% Set filetype and sampling rate we will use.
filetype = SA.filetype;
default_Fs = SA.samprate;

% Set the paths we will use.
path_songfile = SA.rawdatadir;
path_notefile = SA.resultsdir;
path_filtfile = SA.tempdir;
path_spectfile = SA.tempdir;

% Now for something completely different (bdw) ...
% We make a list of files to process so we can do random access.

if SA.single
  SA.song_list{1} = SA.loadfile;
  SA.n_songs = 1
else
  fid = fopen(SA.loadfile,'rt');
  if fid == -1
    error(['Error opening batch file: ', SA.loadfile])
  end
  ifile = 0
  while 1
    fname = fgetl(fid);
    % Check for whitespace here!
    if (~ischar(fname)) | isempty(fname) | isspace(fname)
      disp('End of batch file reached.')
      break
    end
    ifile = ifile+1;
    SA.song_list{ifile} = fname
  end
  SA.n_songs = ifile 
  fclose(fid)
end
setappdata(h_main, 'SA_Data', SA);

%
%main program: cycle through sound files until end or quit command
%

% This indicates which song we are looking at.
i_song = 1

while 1
  %reset values of onsets, offsets, labels, filtered song, spectrogram
  onsets = [];
  offsets = [];
  labels = [];
  filtsong = [];
  spect = [];
  
  %values for threshold, Fs, min_int, min_dur, sm_win will be maintained
  %unless new values are read in from a note file
  
  %get soundfile name
  if i_song > SA.n_songs
    soundfile = '';
    disp('End of songfiles')
    break
  end
  soundfile = SA.song_list{i_song};

  % We allow for the possibility of explicit relative path names in the
  % batch file
  [path_batchsongfile,soundfilename,ext,ver] = fileparts(soundfile);   
  % This is the sound file name with no path spec.
  soundfile = [soundfilename,ext,ver];

  %if soundfile name ends in '.filt' then strip '.filt' from the end
  % so that batch files which list filtfiles can be read
  filt_idx=findstr('.filt',soundfile);
  if ~isempty(filt_idx)
    soundfile=soundfile(1:filt_idx(length(filt_idx))-1);
  end

  % If it looks like the song file has an explicit path name,
  % we ignore the user's specification of the base song directory
  % and use the batch file specification instead.

  if ~isempty(path_batchsongfile)
    % test for windows drive spec. Hack!
    drivespec = sscanf(path_batchsongfile,'%c',2);
    if strncmp(path_batchsongfile,filesep,1) | strcmp(drivespec(2),':')
      % Batch file path is absolute
      if ~isempty(SA.rawdatadir)
	if ~strcmp(SA.rawdatadir,path_batchsongfile)
	  warning(['Batch file song path overrides user specification: user', ...
		' path = ', SA.rawdatadir, ' batch path: ', path_batchsongfile]);
	end
      end
      path_songfile = path_batchsongfile;
    else
      % Batch file path is relative; append user specified path to it
      path_songfile = [SA.rawdatadir, filesep, path_batchsongfile] 
    end
  end
  if ~exist(path_songfile,'dir')
    disp('Song file directory does not exist. Skipping file.')
    continue
  end
  
  %get note_file name
  note_file=[soundfile,'.not.mat'];  
  %if notefile exists, get it
  if exist(fullfile(path_notefile, note_file),'file')
    disp('loading notefile...');     
    load(fullfile(path_notefile, note_file));
    default_Fs = Fs;   %if  Fs was read, it is the new default value
  end
  
  %if file was loaded then: 
  % Fs, threshold, min_int, min_dur, sm_win, onsets, offsets and labels should be defined by notefile values
  % otherwise onsets, offsets and labels are empty, but other values are set at defaults
  
  
  %if filtered song exists, get it, otherwise read and filter song, save filtered song if flag is set
  %note: if filtsong exists, it is unneccessary to have actual songfile around    

  % If user specifies filt type of data then filt files are not
  % found in the tempdir but the data dir.
  if strcmp(filetype,'filt')
    path_filtfile = path_songfile;
  end
  
  filt_file=fullfile(path_filtfile, [soundfile, '.filt']);
  if exist(filt_file,'file')     
    disp('loading filtered song...');
    [filtsong, Fs] = read_filt(filt_file);
    default_Fs = Fs;   %if Fs was read, it is the new default value (takes precedence over notefile if discrepant)
  else    
    [rawsong,Fs]=evsoundin(path_songfile, soundfile, filetype);
    Fs
    %unless Fs has been read from the file use the default value which was
    %set by the user or previously read in from a note file
    if Fs == -1; Fs = default_Fs; end
    if do_filt
      %filter soundfile
      disp('filtering song...');
      filtsong=bandpass(rawsong,Fs,F_low,F_high,filter_type);
      %save soundfile if flag is set
      if save_filt == 1     
	disp('saving filtered song...')
	write_filt(filt_file, filtsong, Fs);
      end
    else
      filtsong = rawsong;      
    end
    clear rawsong;
  end
  
  %if spectrogram is going to be saved or displayed
  % then see if it already exists, if not, calculate it and save if flag is set
  % Even if only printing the figure we still need the spectrogram! (bdw)
  if (batch_disp == 1) | (save_spect == 1) | (batch_print)
    spect_file=fullfile(path_spectfile, [soundfile,'.spect']);
    if exist(spect_file,'file')     
      disp('reading spectrogram...')
      [idx_spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max] = read_spect(spect_file);
      time_spect = [t_min, t_max];
      freq_spect = [f_min, f_max];     
    else
      %calculate spectrogram
      disp('calculating spectrogram...');
      %first calculate nfft and noverlap
      nfft=round(Fs*spect_win_dur/1000);
      nfft = 2^nextpow2(nfft);
      spect_win = hanning(nfft);
      noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points       
      %now calculate spectrogram
      [spect, freq, time] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
      idx_spect=scale_spect(spect);  %calculate index array for spectrogram
      f_min = freq(1);
      f_max = freq(length(freq));
      freq_spect = [f_min, f_max];
      t_min = time(1)*1000; %convert to ms
      t_max = time(length(time))*1000; %convert to ms
      %adjust time axis for spectrogram offset (1/2 window duration in ms)
      t_min = t_min + 0.5*nfft*1000/Fs;  
      t_max = t_max + 0.5*nfft*1000/Fs;  
      time_spect = [t_min, t_max];                
      %save spectrogram if flag is set        
      if save_spect == 1     
	disp('saving spectrogram...')
	write_spect(spect_file, idx_spect, nfft, spect_win, noverlap, t_min, t_max, f_min, f_max);
      end
    end      
  end   
  
  disp('calculating power and smoothing...')
  %calculate square of signal: proportional to power
  squared_song = filtsong.^2;
  
  %smooth the rectified song
  len=round(Fs*sm_win/1000);                      
  h=ones(1,len)/len;
  smooth=conv(h, squared_song);
      save tmp.mat smooth
  offset=round((length(smooth)-length(filtsong))/2); %get rid of convolution induced offset
  smooth=smooth(1+offset:length(filtsong)+offset);
  
  % Zero crossings
  if do_zc
    zcwinlen=round(Fs*zcwin/1000);                      
    % This gets average # of zero crossings per sample
    zcwindow=ones(1,zcwinlen)/zcwinlen;
    zc = zero_crossings(filtsong,zcwindow);
    offsetzc=round((length(zc)-length(filtsong))/2); %get rid of convolution induced offset
    % convert to a frequency estimate
    zc = (Fs/2)*zc(1+offsetzc:length(filtsong)+offsetzc); 
  end
  
  %if threshold is undefined, calculate threshold and segment: will only be true for first file
  if threshold == 0
    disp('calculating default threshold...')
    %get threshold value
    len1=round(Fs*10/1000);
    [m, i]=max(smooth(len1:length(smooth)-len1));
    peak=sum(smooth(i-len1:i+len1))/(2*len1+1);
    threshold=round(peak)/1000;
  end   
  
  
  %if values were not loaded for onsets offsets and labels, calculate them, and save them if flag is set
  if isempty(onsets)   
    disp('segmenting song...')
    [onsets, offsets] = segment(smooth, Fs, min_int, min_dur, threshold);
    %set label values
    labels= num2str(zeros(size(onsets)))';
    %save note data if flag is set
    if save_notes == 1; 
      save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win);
    end 
  end  

  %fix for matlab 5.x compatability???
  labels=makerow(labels);
  onsets=makecol(onsets);
  offsets=makecol(offsets);
  
  %display only if flag is set
  if batch_disp == 1 | batch_print == 1
    
    %resample song for more rapid display
    disp('resampling for display...');
    step=round(Fs/d_Fs);
    d_Fs=Fs/step;
    d_song=smooth(1:step:length(smooth));
    if do_zc
      d_zc=zc(1:step:length(zc));
    end
    
    %set initial values of display variables  
    time_d_song=[0:length(d_song)-1]*1000/d_Fs;  %vector for time axis
    g_xmin=time_d_song(1);                                    %initially display everything
    g_xmax=time_d_song(length(time_d_song));
    g_ymin=0;
    g_ymax=max(d_song)*1.2;
    % The following wierdness is a workaround because of matlab bug for image display
    %see disp_spect for a bit more info
    %basic idea is to invert the spectrogram before displaying and to label y axis with negative 
    %frequency values
    %   g_fmin=(-1)*F_high;
    %   g_fmax=0;  
    g_fmin=0;
    g_fmax=F_high;  
    
    % Turn on figure if in batch or interact modes
    if batch_disp
      %display song
      disp('displaying song...')

      %%%%%%%%%%%%%%%%%%%first time we see the main window%%%%%%%%%%%%%%%
      % Turn on main window
      set(h_main,'Visible','on')
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %could retain old figure here
    %if isempty(h_main);
      %% In an axis tagged "SAAxis", set the pointer to a cross
      %h_main=figure('Tag','SAMain','WindowButtonMotionFcn','centercb(''center_pointer'')','WindowButtonDownFcn','centercb(''center'')');
    %end 
    %see if this helps w/ gradual slowing of program as files are processed
    %works but doesn't preserve window size
    %delete(h_main);
    %h_main=figure;
    
    % Set up plot positions
    pos_axes_left = 0.1300
    pos_amp_axes_bot = 0.1100
    width_axes = 0.7750
    height_amp_axes = 0.3439
    pos_spect_axes_bot = 0.5811;
    height_spect_axes = height_amp_axes;   
    height_label_axes = 0.2*height_amp_axes;
    amptop = pos_amp_axes_bot + height_amp_axes;
    pos_label_axes_bot = amptop + (pos_spect_axes_bot - amptop - height_label_axes)/2;
    
    % h_main_amp=subplot(2,1,2);
    h_main_amp=subplot('Position', [pos_axes_left pos_amp_axes_bot width_axes height_amp_axes]);
    % Tag for identifying main axes 
    set(h_main_amp,'Tag','SAAxis')
    %store the global values of data limits with axes
    set(h_main_amp,'userdata',[g_xmin g_xmax g_ymin g_ymax]);
    if xrange_start == -1 | g_xmin + xrange_start > g_xmax
      %set the current display to show all the data
      set(h_main_amp,'xlim',[g_xmin g_xmax]);
    else
      set(h_main_amp,'xlim',[g_xmin g_xmin + xrange_start]);
    end
    set(h_main_amp,'ylim',[g_ymin g_ymax]);
    
    % Here make the subplot for labels
    h_main_labels=subplot('Position', [pos_axes_left pos_label_axes_bot width_axes height_label_axes]);
    % Tag for identifying main axes 
    set(h_main_labels,'Tag','SAAxis','Visible','off')
    %store the global values of data limits with axes
    set(h_main_labels,'userdata',[g_xmin g_xmax 0 1]);
    if xrange_start == -1 | g_xmin + xrange_start > g_xmax
      %set the current display to show all the data
      set(h_main_labels,'xlim',[g_xmin g_xmax]);
      left = g_xmin;
      right = g_xmax;
    else
      set(h_main_labels,'xlim',[g_xmin g_xmin + xrange_start]);
      left = g_xmin;
      right = g_xmin + xrange_start;
    end
    set(h_main_labels,'ylim',[0 1]);
    h_patch = patch('xdata',[left left right right],'ydata',[0.25 0.75 ...
	  0.75 0.25], 'FaceColor', 'w', 'LineStyle', 'none');
    
    [h_amp_plot, h_labels] = disp_song;

    %display spectrogram
    if batch_disp
      disp('displaying spectrogram...')
    end
    % h_main_spect=subplot(2,1,1);
    h_main_spect=subplot('Position',[pos_axes_left pos_spect_axes_bot width_axes height_spect_axes]);
    % Tag for identifying main axes 
    set(h_main_spect,'Tag','SAAxis')
    %   set(h_main_spect,'ydir', 'reverse');
    if xrange_start == -1 | g_xmin + xrange_start > g_xmax
      %initially display all data
      set(h_main_spect,'xlim',[g_xmin g_xmax]);
    else
      set(h_main_spect,'xlim',[g_xmin g_xmin + xrange_start]);
    end
    set(h_main_spect,'ylim',[g_fmin g_fmax]);
    %store max and min data values in user data
    set(h_main_spect,'userdata',[g_xmin g_xmax g_fmin g_fmax]);
    title(soundfile);
    disp_idx_spect(idx_spect, time_spect, freq_spect, spect_floor, ...
	spect_ceil, spect_range, spect_cmap_name, spect_gamma_type);

    if do_zc
      hold on
      h_zc_plot = plot(time_d_song, d_zc, 'g-');
      set(h_zc_plot,'Tag','ZCPlot');
      hold off
    end
    
    subplot(h_main_amp);
    
    %set some printing defaults
    set_print_vals(h_main, pr_orient, pr_papertype, pr_x_per_paperinch, g_ymax-g_ymin); %y value is superceded on next line
    temp = get(h_main, 'paperposition');
    temp(4) = pr_y_paperinches;
    set(h_main, 'paperposition', temp);
    
    %reset data status for any open windows
    windows = get(0,'children');
    for i = 1:length(windows)
      userdata = get(windows(i),'userdata');
      if length(userdata > 0)
	if userdata(imain_win_code) == main_win_code;
	  userdata(idata_status) = 0;
	  set(windows(i),'userdata',userdata);
	end
      end
    end
    
    %store file and pathnames as hidden text
    %store these with spectrogram axes for now: note this will be a problem if spectrogram axis
    %is cleared
    subplot(h_main_spect);
    h_soundfile = text(-100, -100, soundfile, 'visible', 'off');
    h_path_songfile = text(-100, -100, path_songfile, 'visible', 'off');                       
    h_path_notefile = text(-100, -100, path_notefile, 'visible', 'off');   
    h_path_filtfile = text(-100, -100, path_filtfile, 'visible', 'off'); 
    h_path_spectfile = text(-100, -100, path_spectfile, 'visible', 'off'); 
    
    %set userdata values
    userdata(imain_win_code) = main_win_code;
    userdata(idata_status) = 3;  %see make_current for meaning of code
    userdata(ih_main_amp) = h_main_amp;
    userdata(ih_main_spect) = h_main_spect;
    userdata(ih_fname) = h_soundfile; 
    userdata(ih_path_songfile) = h_path_songfile;  
    userdata(ih_path_notefile) = h_path_notefile;
    userdata(ih_path_filtfile) = h_path_filtfile;
    userdata(ih_path_spectfile) = h_path_spectfile;
    userdata(ih_amp_plot) = h_amp_plot;    
    userdata(ih_main_labels) = h_main_labels;
    if do_zc
      userdata(ih_zc_plot) = h_zc_plot;
    end  
    subplot(h_main_amp);   
    
    %store info in userdata for retrieval of handles etc
    set(h_main,'userdata',userdata);
    
    %build ui controls: sets up all the interface controls in current figure
    uisongcontrols('initialize');
    uispectcontrols('initialize');  
    
    %if in batch mode and want to display, need to pause here
    if interact == 0 & batch_disp == 1
      pause(.3);
    end
    
  end   %end of building window if in display mode
  
  %if batch_print is set, print window
  if batch_print == 1
    %set size of current window to fit page
    set(gcf,'paperorientation', pr_orient);
    %save unit values and set to inches
    t_paperunits = get(gcf,'PaperUnits');
    set(gcf,'PaperUnits','inches');
    t_fig_units = get(gcf,'Units');
    set(gcf,'Units','inches');
    %set figure on screen to match page size
    % first move the figure to bottom left to avoid going off
    % screen as well as crash under linux that happens if it does!
    movegui(gcf,'southwest');
    position = get(gcf,'Position');
    position(3)=11;
    position(4)= pr_y_paperinches;
    set(gcf,'Position',position);    

    %set horizontal scale of current window to desired value
    set_scale_vals(gcf, pr_x_per_paperinch, pr_y_paperinches);
%    movegui(gcf,'onscreen')

    %set paperposition values
    position(1) = pr_left_margin;
    position(2) = pr_bottom_margin;
    set(gcf,'paperposition',position);
    %print & shift until done with spectrogram
    more = 1;
    while more == 1
      print -dps2 -noui
      subplot(h_main_amp);
      more = move_right(win_percent);
      subplot(h_main_labels);
      move_right(win_percent);
      subplot(h_main_spect);
      move_right(win_percent);
      subplot(h_main_amp);
    end
    %reset units  
    set(gcf,'PaperUnits',t_paperunits);
    set(gcf,'Units',t_fig_units);
  end
    
  %sit idle until need to build another window or in batch mode
  while get_next==0 & interact == 1
    pause(.5)
  end
  
  %reset
  get_next = 0;

  % Got command to go to the next song.
  i_song = i_song+1;
  
end

if ~batch_disp
  % if figure is invisible, close it at end
  disp('Closing figure...');
  close(gcf);
end
% remove the temporary sub path(?)
rmpath(subpath);

disp('Done!');

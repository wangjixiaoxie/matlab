function varargout = filecb(action,varargin)

disp(action);

switch action
  case 'loadfile'
    if nargin == 1
      h_main = gcbf;
    elseif nargin == 2
      h_main = varargin{1};
    end

    % Start up the file loading dialog    
    fig = uisasetup;
    % Set up defaults
    SALoad = getappdata(fig,'SALoad_Data');
    SALoad.handles.SALoad = fig;
    SALoad.handles.SAMain = h_main;
    SA = getappdata(h_main,'SA_Data')
    setappdata(fig,'SALoad_Data',SALoad)

    % Set some defaults, depending on whether a file is already loaded.
    if SA.file_loaded == 0
      filecb('filetype','init',fig);
      SALoad = getappdata(fig,'SALoad_Data');
      SALoad.samprate = SA.samprate;
      set(SALoad.handles.SALoadSamplingRateEd,'String',num2str(SALoad.samprate))
      set(SALoad.handles.SALoadFileCkb,'Value',0);
      set(SALoad.handles.SALoadRawdirEd,'Enable','on')
      set(SALoad.handles.SALoadResultsdirEd,'String',pwd)
      set(SALoad.handles.SALoadTempdirEd,'String',pwd)
    else
      filecb('filetype','set',fig);
      SALoad = getappdata(fig,'SALoad_Data');
      set(SALoad.handles.SALoadSamplingRateEd,'String',num2str(SALoad.samprate))
      set(SALoad.handles.SALoadFileCkb,'Value',SA.single);
      if SA.single == 1
	set(SALoad.handles.SALoadRawdirEd,'Enable','off','BackgroundColor',[0.8,0.8,0.8])      
      end
      set(SALoad.handles.SALoadFileEd,'String',SA.loadfile);
      set(SALoad.handles.SALoadRawdirEd,'String',SA.rawdatadir)
      set(SALoad.handles.SALoadResultsdirEd,'String',SA.resultsdir)
      set(SALoad.handles.SALoadTempdirEd,'String',SA.tempdir)
    end
    platf = computer;
    % Turn off java directory browser on windows until matlab 6.5 upgrade
    if strcmp(platf,'PCWIN') & str2num(version('-release')) <= 12.1 
      set(SALoad.handles.SALoadBrowseRawdirPb,'Enable','off');
      set(SALoad.handles.SALoadBrowseResultsdirPb,'Enable','off');
      set(SALoad.handles.SALoadBrowseTempdirPb,'Enable','off');
    end
    setappdata(fig,'SALoad_Data',SALoad)

    varargout{1} = fig;
    % Process callbacks until dialog is closed.
    uiwait(fig);
    
  case 'single'
    SALoad = getappdata(gcbf,'SALoad_Data');
    singleflag = get(SALoad.handles.SALoadFileCkb,'Value');
    if singleflag == 1
      set(SALoad.handles.SALoadRawdirEd,'Enable','off','BackgroundColor',[0.8,0.8,0.8])
    else
      set(SALoad.handles.SALoadRawdirEd,'Enable','on','BackgroundColor',[1.0,1.0,1.0])
    end
  
  case 'filetype'
    if nargin == 3
      fig = varargin{2};
      myaction = varargin{1}
    elseif nargin == 2
      fig = varargin{1};
    elseif nargin == 1
      fig = gcbf;
    end
    SALoad = getappdata(fig,'SALoad_Data');
    
    if exist('myaction','var') & strcmp(myaction,'init')    
      % set defaults
      SALoad.filetype = 'w';
      SALoad.filetype_strs = {'w', 'obs0r', 'obs1r', 'd', 'b', 'filt', 'foo'};
      SALoad.filetype_tips = {'WAV file', 'Observer (last channel)', ...
	    'Observer (next to last channel)', 'dcp file', ...
	    'binary (mac) format', 'filt file', 'foosong file'}
      filetypeidx = find(strcmpi(SALoad.filetype_strs,SALoad.filetype));
      if ~isempty(filetypeidx)
	set(SALoad.handles.SALoadFileTypePu,'Value',filetypeidx,'TooltipString', SALoad.filetype_tips{filetypeidx});
      end
      setappdata(fig,'SALoad_Data',SALoad)

    elseif exist('myaction','var') & strcmp(myaction,'set')    
      filetypeidx = find(strcmpi(SALoad.filetype_strs,SALoad.filetype));
      if ~isempty(filetypeidx)
	set(SALoad.handles.SALoadFileTypePu,'Value',filetypeidx,'TooltipString', SALoad.filetype_tips{filetypeidx});
      end
    elseif exist('myaction','var') & strcmp(myaction,'get')    
      filetypeidx = get(SALoad.handles.SALoadFileTypePu,'Value');
      SALoad.filetype = SALoad.filetype_strs{filetypeidx}
      setappdata(fig,'SALoad_Data',SALoad)    
    end    
    
    filetypeidx = get(SALoad.handles.SALoadFileTypePu,'Value');
    set(SALoad.handles.SALoadFileTypePu, 'TooltipString', SALoad.filetype_tips{filetypeidx});
    filetype = SALoad.filetype_strs{filetypeidx};
    if strncmp(filetype,'w',1) | strncmp(filetype,'filt',4)
      set(SALoad.handles.SALoadSamplingRateEd,'Enable','off','BackgroundColor',[0.8,0.8,0.8])
    else
      set(SALoad.handles.SALoadSamplingRateEd,'Enable','on','BackgroundColor',[1.0,1.0,1.0])
    end
    
  case 'browseload'
    % Select a load file and set the associated edit box
    SALoad = getappdata(gcbf,'SALoad_Data');
    dlgtitle = 'Select file';
    loadfile = filecb('browse',dlgtitle)
    if ~isempty(loadfile)
      varargout{1} = loadfile;
      set(SALoad.handles.SALoadFileEd,'String',loadfile);
      [pname,fname,ext] = fileparts(loadfile);
      if exist(pname,'dir')
	set(SALoad.handles.SALoadRawdirEd,'String',pname);
      end	
    end
  
  case 'browserawdir'
    % Select a directory and set the associated edit box
    dlgtitle = 'Select raw data directory';
    SALoad = getappdata(gcbf,'SALoad_Data');
    h_edit = SALoad.handles.SALoadRawdirEd;
    
    initdir = get(h_edit,'String');
    if isempty(initdir)
      initdir = pwd;
    end
    choosedir = filecb('browsedir',initdir,dlgtitle)

    if ~isempty(choosedir)
      set(h_edit,'String', choosedir);
      varargout{1} = choosedir;
    end

  case 'browseresultsdir'
    % Select a directory and set the associated edit box
    dlgtitle = 'Select directory for results';
    SALoad = getappdata(gcbf,'SALoad_Data');
    h_edit = SALoad.handles.SALoadResultsdirEd; 

    initdir = get(h_edit,'String');
    if isempty(initdir) | ~exist(initdir,'dir')
      initdir = pwd;
    end
    choosedir = filecb('browsedir',initdir,dlgtitle)

    if ~isempty(choosedir)
      set(h_edit,'String', choosedir);
      varargout{1} = choosedir;
    end

  case 'browsetempdir'
    % Select a directory and set the associated edit box
    dlgtitle = 'Select directory for temp files';
    SALoad = getappdata(gcbf,'SALoad_Data');
    h_edit = SALoad.handles.SALoadTempdirEd;

    initdir = get(h_edit,'String');
    if isempty(initdir) | ~exist(initdir,'dir')
      initdir = pwd;
    end
    choosedir = filecb('browsedir',initdir,dlgtitle)

    if ~isempty(choosedir)
      set(h_edit,'String', choosedir);
      varargout{1} = choosedir;
    end

  case 'browse'
    % Generic file selector
    if nargin == 2 & ischar(varargin{1})
      dlgtitle = varargin{1};
      [fname, pname] = uigetfile('*',dlgtitle)
      if ~isequal(fname,0) & ~isequal(pname,0)
	varargout{1} = fullfile(pname,fname)	
      end
    else 
      error('Incorrect argument list to filecb.')
    end
      
  case 'browsemat'
    % Generic mat file selector
    if nargin == 2 & ischar(varargin{1})
      dlgtitle = varargin{1};
      [fname, pname] = uigetfile('*.mat',dlgtitle)
      if ~isequal(fname,0) & ~isequal(pname,0)
	varargout{1} = fullfile(pname,fname)	
      end
    else 
      error('Incorrect argument list to filecb.')
    end

  case 'putmat'
    % Generic save mat file selector
    if nargin == 2 & ischar(varargin{1})
      dlgtitle = varargin{1};
      [fname, pname] = uiputfile('*.mat',dlgtitle)
      if ~isequal(fname,0) & ~isequal(pname,0)
	varargout{1} = fullfile(pname,fname)	
      end
    else 
      error('Incorrect argument list to filecb.')
    end

  case 'browsedir'
    % Generic directory selector
    if nargin == 1
      if str2num(version('-release')) > 12.1 
	pname = uigetdir 
      else
	pname = uigetdir_java 
      end
    elseif nargin == 2 & ischar(varargin{1})
      initdir = varargin{1};
      if str2num(version('-release')) > 12.1 
	pname = uigetdir(initdir) 
      else
	pname = uigetdir_java(initdir)
      end
    elseif nargin == 3 & ischar(varargin{1}) & ischar(varargin{2})
      initdir = varargin{1};
      dlgtitle = varargin{2};
      if str2num(version('-release')) > 12.1 
	pname = uigetdir(initdir,dlgtitle)
      else
	pname = uigetdir_java(initdir,dlgtitle)
      end
    else 
      error('Incorrect argument list to filecb.')
    end
    varargout{1} = pname;
    
  case 'ok'
    SALoad = getappdata(gcbf,'SALoad_Data');
    filecb('apply')
    uiresume(gcbf);
    close(SALoad.handles.SALoad); 
    
  case 'apply'

    filecb('filetype','get',gcbf);
    SALoad = getappdata(gcbf,'SALoad_Data');
    SA = getappdata(SALoad.handles.SAMain,'SA_Data');
    SA.single = get(SALoad.handles.SALoadFileCkb,'Value');

    SA.filetype = SALoad.filetype;
    if ~(strncmp(SA.filetype,'w',1) | strncmp(SA.filetype,'filt',4))
      SA.samprate = str2num(get(SALoad.handles.SALoadSamplingRateEd,'String'));
    end
    SA.resultsdir = get(SALoad.handles.SALoadResultsdirEd,'String')
    SA.tempdir = get(SALoad.handles.SALoadTempdirEd,'String')
    
    loadfile = get(SALoad.handles.SALoadFileEd,'String');
    loadfail = 0;
    if ~isempty(loadfile) & exist(loadfile,'file')
      SA.loadfile = loadfile
      SA.file_loaded = 1;
    else
      loadfail = 1;
      filecb('error_file',loadfile)
    end

    if SA.single == 0
      adir = get(SALoad.handles.SALoadRawdirEd,'String')
      if exist(adir,'dir')
	SA.rawdatadir = adir
      else
	filecb('error_dir',dir)	
      end
    end
    adir = get(SALoad.handles.SALoadResultsdirEd,'String')
    if exist(adir,'dir')
      SA.resultsdir = adir
    else
      filecb('error_dir',adir)
    end
    adir = get(SALoad.handles.SALoadTempdirEd,'String')
    if exist(adir,'dir')
      SA.tempdir = adir
    else
      filecb('error_dir',adir)
    end
    setappdata(SALoad.handles.SAMain,'SA_Data',SA);    
      
  case 'cancel'
    SALoad = getappdata(gcbf,'SALoad_Data');
    uiresume(gcbf);
    close(SALoad.handles.SALoad); 

  case 'error_dir'
    if nargin == 2 & ischar(varargin{1})
      warnstr = ['Directory ', varargin{1}, ' does not exist. Nothing Loaded.'];
    else
      warnstr = 'Directory does not exist. Nothing Loaded.';
    end
    warndlg(warnstr,'Warning')
    
  case 'error_file'
    if nargin == 2 & ischar(varargin{1})
      warnstr = ['File ', varargin{1}, ' does not exist. Nothing Loaded.'];
    else
      warnstr = 'File does not exist. Nothing Loaded.';
    end
    warndlg(warnstr,'Warning')
    
 otherwise
   error('Unknown file callback action.')
   
end
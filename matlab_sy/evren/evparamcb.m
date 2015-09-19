function varargout = evparamcb(action,varargin)

disp(action);

switch action
  case 'construct'
    if nargin == 1
      h_main = gcbf;
    elseif nargin >= 2
      h_main = varargin{1};
    end
    if nargin == 3
      config_file = varargin{2};
    end
    
    % This will be the tabbed dialog
    % tabstrings = {' General ',' Segmentation/Notes ',' Filtering',' Spectrogram ',' Printing ', ' Misc '};

    % tabdims = tabdlg('tabdims', tabstrings);
    % [fig, sheetpos] = tabdlg('create', tabstrings, tabdims, ...
    %	'tabcbk', [600, 425], [25,10,25,100],1) 

    fig = uisasetparams;
    SASet = getappdata(fig,'SASet_Data');
    SASet.handles.SASet = fig;
    SASet.handles.SAMain = h_main;
    SA = getappdata(SASet.handles.SAMain,'SA_Data');

    SASet.tabstrings = {' General ',' Segmentation/Notes ',' Filtering ',' Spectrogram ',' Printing ', ' Misc '};
    ntabs = length(SASet.tabstrings);
    
    % Note all tags on a given tabbed page must contain the associated
    % tabfield as a substring!
    SASet.tabfields = {'General','Segment','Filter','Specgram','Print', 'Misc'};    
    % Initialize structure fields
    for ifld = 1:ntabs
      eval(['SASet.', SASet.tabfields{ifld}, '.handles = [];']);
    end
    
    % Load in all the tags and segregate handles according to tab page
    all_tags = get(get(fig,'Children'),'Tag');    
    for it = 1:length(all_tags)
      ifld = 1;
      while ifld <= length(SASet.tabfields) & isempty(findstr(all_tags{it}, SASet.tabfields{ifld}))
	ifld = ifld + 1;
      end
      if ifld > ntabs
	% disp(['Tag: ', all_tags{it}, ' not associated with a tab page.'])
	continue;
      else
	% disp(all_tags{it})
	eval(['SASet.', SASet.tabfields{ifld}, '.handles = [SASet.', ...
	      SASet.tabfields{ifld}, '.handles, SASet.handles.', ...
	      all_tags{it}, '];']);
      end
    end

    % Make all pages invisible initially.
    for ifld = 1:ntabs
      eval(['set(SASet.', SASet.tabfields{ifld}, '.handles, ''Visible'',''off'');']);
    end
    
    setappdata(fig,'SASet_Data',SASet);

    % Create tabs and reset the tab dialog to default page.
    default_page = 1;
    dummy_page = ntabs;
    evparamcb('create_tabs', fig, dummy_page);
    tabdlg('tabpress',fig, SASet.tabstrings{default_page}, default_page);
    % evparamcb('reset_tabdlg', fig, SASet.tabstrings{default_page}, default_page);

    % Set some defaults, depending on whether parameters have
    % already been set.
    if SA.params_loaded == 0
      evparamcb('use_defaults',fig);
      if exist('config_file','var') & exist(config_file,'file')
	evparamcb('restore_config',config_file,fig);
      end
    else
      % Here set params to current values instead of defaults
      % Current values must be stored in the main UI appdata
      evparamcb('use_current',fig);
    end

    % After set up make the dialog visible.
    set(fig,'Visible','on');
    varargout{1} = fig;
    % Process callbacks until dialog is closed.
    uiwait(fig);
  
  case 'use_defaults'
    % These are the "factory" defaults.
    if nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');
    SA = getappdata(SASet.handles.SAMain,'SA_Data');
    % SA.param_defaults should be a cell array of defaults for all user
    % parameters, stored in the main UI.  
    
    % Set all parameters to defaults
    params = fieldnames(SA.param_defaults);   
    for ip=1:length(params)
      eval(['SASet.params.', params{ip}, ' = SA.param_defaults.', params{ip},';']);
    end
    setappdata(fig,'SASet_Data',SASet);

    % Actual control values set here.
    evparamcb('initcontrols',fig);

  case 'use_current'
    % These are the current values stored in the main UI.
    if nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');
    SA = getappdata(SASet.handles.SAMain,'SA_Data');
    % SA.params should be a cell array of current values for all user
    % parameters, stored in the main UI.  
    
    % Set all parameters to defaults
    params = fieldnames(SA.param_defaults);   
    for ip=1:length(params)
      eval(['SASet.params.', params{ip}, ' = SA.params.', params{ip},';']);
    end
    setappdata(fig,'SASet_Data',SASet);

    % Actual control values set here.
    evparamcb('initcontrols',fig);

  case 'initcontrols'
    fig = varargin{1};
    % Reinitalize all controls
    evparamcb('general_params','init', fig);
    evparamcb('segment_params','init', fig);
    evparamcb('filter_params','init', fig);
    evparamcb('specgram_params','init', fig);
    evparamcb('print_params','init', fig);
    evparamcb('misc_params','init', fig);

  case 'get_params'
    fig = varargin{1};
    % Update all appdata to current settings in dialog
    evparamcb('general_params','get', fig);
    evparamcb('segment_params','get', fig);
    evparamcb('filter_params','get', fig);
    evparamcb('specgram_params','get', fig);
    evparamcb('print_params','get', fig);
    evparamcb('misc_params','get', fig);

    SASet = getappdata(fig,'SASet_Data');

    % Update current params list
    ntabs = length(SASet.tabstrings);
    SA = getappdata(SASet.handles.SAMain,'SA_Data');
    allparams = fieldnames(SA.param_defaults);
    for ifld = 1:ntabs
      eval(['tabparams = fieldnames(SASet.', ...
	    SASet.tabfields{ifld}, ');']);
      ntabparams = length(tabparams);
      for ip = 1:ntabparams
	if any(strcmp(tabparams{ip},allparams))
	  % disp(['SASet.params.', tabparams{ip}, ' = SASet.', ...
	  %	SASet.tabfields{ifld}, '.', tabparams{ip}])	  
	  eval(['SASet.params.', tabparams{ip}, ' = SASet.', ...
		SASet.tabfields{ifld}, '.', tabparams{ip},';']);	  
	end
      end
    end
    setappdata(fig,'SASet_Data',SASet);      

  case 'restore_config'
    % Restore configuration from a .mat file
    if nargin == 3
      fig = varargin{2};
      config_file = varargin{1};
    elseif nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');
    SA = getappdata(SASet.handles.SAMain,'SA_Data');

    if ~exist('config_file','var') | ~exist(config_file,'file')
      dlgtitle = 'Select configuration file';
      config_file = filecb('browsemat',dlgtitle)
    end
    % Empty file means user cancelled the restore?
    if ~isempty(config_file)
      [cfgpathstr,cfgname,cfgext,cfgvers] = fileparts(config_file); 
      if strcmpi(cfgext,'.mat')
	myconfig = load(config_file);
	if isfield(myconfig,'SASet')
	  if ~isfield(myconfig.SASet,'params')
	    warn(['Invalid config file: ', config_file])
	    return
	  end
	else
	  warn(['Invalid config file: ', config_file])
	  return
	end
	% File passes all checks
	SASet.config_file = config_file;

	myparams = fieldnames(myconfig.SASet.params);
	allparams = fieldnames(SA.param_defaults);
	for ip=1:length(myparams)
	  % disp(['myparam: ', myparams{ip}])
	  if any(strcmp(myparams{ip},allparams)) 
	    % disp(['SASet.params.', myparams{ip}, ' = myconfig.SASet.params.', myparams{ip}]);
	    eval(['SASet.params.', myparams{ip}, ' = myconfig.SASet.params.', myparams{ip},';']);
	  end
	end
	setappdata(fig,'SASet_Data',SASet);      

	% Update control values.
	evparamcb('initcontrols',fig);
      else
	warn(['Invalid config file: ', config_file])
      end    
    end
    
  case 'save_config'
    % Save config to a .mat file
    if nargin == 3
      fig = varargin{1};
      config_file = varargin{2};
    elseif nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end

    if ~exist('config_file','var')
      dlgtitle = 'Save configuration as';
      config_file = filecb('putmat',dlgtitle)
    end
    % Empty file means user cancelled the save?
    if isempty(config_file)
      return
    end

    % Update SASet structure based on current values in all controls
    evparamcb('get_params',fig);

    SASet = getappdata(fig,'SASet_Data');
    SASet.config_file = config_file;
    setappdata(fig,'SASet_Data',SASet);      

    % OK, now save data into structure
    % We stuff in the whole SASet structure.
    save(config_file,'SASet');
        
  case 'ok'
    if nargin == 2      
      fig = varargin{1};
    else
      fig = gcbf;
    end
    evparamcb('apply',fig);
    % evparamcb('save_config',fig);
    SASet = getappdata(fig,'SASet_Data');
    uiresume(fig);
    close(SASet.handles.SASet);
  
  case 'apply'
    % 'Apply' means get all current data from controls and
    % set the parameters in the main figure. It could also
    % ask the main UI to update its display based on the new values.
    if nargin == 2      
      fig = varargin{1};
    else
      fig = gcbf;
    end

    % Update SASet structure based on current values in all controls
    evparamcb('get_params',fig);
    
    SASet = getappdata(fig,'SASet_Data');
    SA = getappdata(SASet.handles.SAMain,'SA_Data');
    ntabs = length(SASet.tabstrings);
    allparams = fieldnames(SA.param_defaults);
    for ifld = 1:ntabs
      eval(['tabparams = fieldnames(SASet.', ...
	    SASet.tabfields{ifld}, ');']);
      ntabparams = length(tabparams);
      for ip = 1:ntabparams
	if any(strcmp(tabparams{ip},allparams))
	  % disp(['Matching parameter: ', tabparams{ip}])
	  % disp(['SA.params.', tabparams{ip}, ' = SASet.', ...
	  %	SASet.tabfields{ifld}, '.', tabparams{ip}]);
	  eval(['SA.params.', tabparams{ip}, ' = SASet.', ...
		SASet.tabfields{ifld}, '.', tabparams{ip},';']);
	end
      end
    end
    SA.params_loaded = 1;
    if isfield(SASet,'config_file')
      SA.config_file = SASet.config_file;
    end
    setappdata(SASet.handles.SAMain,'SA_Data', SA);
    
  case 'cancel'
    % Cancel means give up on all changes made to the dialog.
    % Main UI will query the SA.params_loaded variable to see if it
    % has to use the defaults. 
    SASet = getappdata(gcbf,'SASet_Data');
    uiresume(gcbf);
    close(SASet.handles.SASet);
    
  % General parameters
  case 'general_params'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'init')
	SASet.General.save_notes = SASet.params.save_notes;
	SASet.General.save_filt = SASet.params.save_filt;
	SASet.General.save_spect = SASet.params.save_spect;
	SASet.General.interact = SASet.params.interact;
	SASet.General.batch_disp = SASet.params.batch_disp;
	SASet.General.batch_print = SASet.params.batch_print;
	set(SASet.handles.SASetGeneralSaveNotesCkb,'Value',SASet.General.save_notes);
	set(SASet.handles.SASetGeneralSaveFiltCkb,'Value',SASet.General.save_filt);
	set(SASet.handles.SASetGeneralSaveSpectCkb,'Value',SASet.General.save_spect);
	set(SASet.handles.SASetGeneralBatchModeCkb,'Value',~SASet.General.interact);
	set(SASet.handles.SASetGeneralBatchDisplayCkb,'Value',SASet.General.batch_disp);
	set(SASet.handles.SASetGeneralBatchPrintCkb,'Value',SASet.General.batch_print);
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
	evparamcb('batch_mode','set',fig);
      
      elseif strcmp(varargin{1},'set')
	set(SASet.handles.SASetGeneralSaveNotesCkb,'Value',SASet.General.save_notes);
	set(SASet.handles.SASetGeneralSaveFiltCkb,'Value',SASet.General.save_filt);
	set(SASet.handles.SASetGeneralSaveSpectCkb,'Value',SASet.General.save_spect);
	set(SASet.handles.SASetGeneralBatchModeCkb,'Value',~SASet.General.interact);
	set(SASet.handles.SASetGeneralBatchDisplayCkb,'Value',SASet.General.batch_disp);
	set(SASet.handles.SASetGeneralBatchPrintCkb,'Value',SASet.General.batch_print);
	evparamcb('batch_mode','set',fig);
      
      elseif strcmp(varargin{1},'get')
	SASet.General.save_notes = get(SASet.handles.SASetGeneralSaveNotesCkb,'Value');
	SASet.General.save_filt = get(SASet.handles.SASetGeneralSaveFiltCkb,'Value');
	SASet.General.save_spect = get(SASet.handles.SASetGeneralSaveSpectCkb,'Value');
        batch_mode = get(SASet.handles.SASetGeneralBatchModeCkb,'Value');
	SASet.General.interact = ~batch_mode;
	SASet.General.batch_disp = get(SASet.handles.SASetGeneralBatchDisplayCkb,'Value');
	SASet.General.batch_print = get(SASet.handles.SASetGeneralBatchPrintCkb,'Value');
	% The 'get' method updates the appdata.
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      end
    end

  case 'batch_mode'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'set')
	batch_mode = ~SASet.General.interact;
      end
    else
      SASet = getappdata(gcbf,'SASet_Data');
      batch_mode = get(SASet.handles.SASetGeneralBatchModeCkb,'Value');
    end
  
    if batch_mode == 0
      set(SASet.handles.SASetGeneralBatchDisplayCkb,'Enable','off');
      set(SASet.handles.SASetGeneralBatchPrintCkb,'Enable','off');
    else
      set(SASet.handles.SASetGeneralBatchDisplayCkb,'Enable','on');
      set(SASet.handles.SASetGeneralBatchPrintCkb,'Enable','on');
    end  
    
  % Segmentation/Notes parameters
  case 'segment_params'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'init')
	SASet.Segment.sm_win = SASet.params.sm_win;
	SASet.Segment.min_int = SASet.params.min_int;
	SASet.Segment.min_dur = SASet.params.min_dur;
	SASet.Segment.max_dur = SASet.params.max_dur;
	SASet.Segment.threshold = SASet.params.threshold;
	set(SASet.handles.SASetSegmentSmWinEd,'String',num2str(SASet.Segment.sm_win));
	set(SASet.handles.SASetSegmentMinIntEd,'String',num2str(SASet.Segment.min_int));
	set(SASet.handles.SASetSegmentMinDurEd,'String',num2str(SASet.Segment.min_dur));
	set(SASet.handles.SASetSegmentThresholdEd,'String',num2str(SASet.Segment.threshold));
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      
      elseif strcmp(varargin{1},'set')
	set(SASet.handles.SASetSegmentSmWinEd,'String',num2str(SASet.Segment.sm_win));
	set(SASet.handles.SASetSegmentMinIntEd,'String',num2str(SASet.Segment.min_int));
	set(SASet.handles.SASetSegmentMinDurEd,'String',num2str(SASet.Segment.min_dur));
	set(SASet.handles.SASetSegmentThresholdEd,'String',num2str(SASet.Segment.threshold));
      
      elseif strcmp(varargin{1},'get')
	SASet.Segment.sm_win = str2num(get(SASet.handles.SASetSegmentSmWinEd,'String'));
	SASet.Segment.min_int = str2num(get(SASet.handles.SASetSegmentMinIntEd,'String'));
	SASet.Segment.min_dur = str2num(get(SASet.handles.SASetSegmentMinDurEd,'String'));
	SASet.Segment.threshold = str2num(get(SASet.handles.SASetSegmentThresholdEd,'String'));
	% The 'get' method updates the appdata.
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      end
    end
    
  % Filter parameters
  case 'filter_params'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'init')
	SASet.Filter.do_filt = SASet.params.do_filt;
	SASet.Filter.F_low = SASet.params.F_low;
	SASet.Filter.F_high = SASet.params.F_high;
	SASet.Filter.filter_type = SASet.params.filter_type;
	set(SASet.handles.SASetFilterDoCkb,'Value',SASet.Filter.do_filt);
	set(SASet.handles.SASetFilterFlowEd,'String',num2str(SASet.Filter.F_low));
	set(SASet.handles.SASetFilterFhighEd,'String',num2str(SASet.Filter.F_high));

	SASet.Filter.filter_type_strs = {'butter', 'hanningfir'};
	filteridx = find(strcmpi(SASet.Filter.filter_type_strs,SASet.Filter.filter_type));
	if ~isempty(filteridx)
	  set(SASet.handles.SASetFilterTypePu,'Value',filteridx);
	end
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
	evparamcb('do_filt','set',fig);
      
      elseif strcmp(varargin{1},'set')
	set(SASet.handles.SASetFilterDoCkb,'Value',SASet.Filter.do_filt);
	set(SASet.handles.SASetFilterFlowEd,'String',num2str(SASet.Filter.F_low));
	set(SASet.handles.SASetFilterFhighEd,'String',num2str(SASet.Filter.F_high));
	filteridx = find(strcmpi(SAset.Filter.filter_type_strs,SASet.Filter.filter_type));
	if ~isempty(filteridx)
	  set(SASet.handles.SASetFilterTypePu,'Value',filteridx);
	end
	evparamcb('do_filt','set',fig);
      
      elseif strcmp(varargin{1},'get')
	SASet.Filter.do_filt = get(SASet.handles.SASetFilterDoCkb,'Value');
	SASet.Filter.F_low = str2num(get(SASet.handles.SASetFilterFlowEd,'String'));
	SASet.Filter.F_high = str2num(get(SASet.handles.SASetFilterFhighEd,'String'));
	SASet.Filter.filter_type = SASet.Filter.filter_type_strs{get(SASet.handles.SASetFilterTypePu,'Value')};
	% The 'get' method updates the appdata.
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      end
    end
    
  case 'do_filt'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'set')
	do_filt = SASet.Filter.do_filt;
      end
    else
      SASet = getappdata(gcbf,'SASet_Data');
      do_filt = get(SASet.handles.SASetFilterDoCkb,'Value');
    end
  
    if do_filt == 0
      set(SASet.handles.SASetFilterTypePu,'Enable','off');
      set(SASet.handles.SASetFilterFlowEd,'Enable','off','BackgroundColor',[0.8,0.8,0.8]);
      set(SASet.handles.SASetFilterFhighEd,'Enable','off','BackgroundColor',[0.8,0.8,0.8]);
    else
      set(SASet.handles.SASetFilterTypePu,'Enable','on');
      set(SASet.handles.SASetFilterFlowEd,'Enable','on','BackgroundColor',[1.0,1.0,1.0]);
      set(SASet.handles.SASetFilterFhighEd,'Enable','on','BackgroundColor',[1.0,1.0,1.0]);
    end  

  % Spectrogram parameters
  case 'specgram_params'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'init')
	SASet.Specgram.spect_win_dur = SASet.params.spect_win_dur;
	SASet.Specgram.spect_overlap = SASet.params.spect_overlap;
	SASet.Specgram.spect_floor = SASet.params.spect_floor;
	SASet.Specgram.spect_ceil = SASet.params.spect_ceil;
	SASet.Specgram.spect_range = SASet.params.spect_range;
	SASet.Specgram.spect_cmap_name = SASet.params.spect_cmap_name;
	SASet.Specgram.spect_gamma_type = SASet.params.spect_gamma_type;
	SASet.Specgram.spect_cmap_name_strs = {'gray', 'hot', 'bone'};
	SASet.Specgram.spect_gamma_type_strs = {'classic', 'brightness'};
	set(SASet.handles.SASetSpecgramWinDurEd,'String',num2str(SASet.Specgram.spect_win_dur));
	set(SASet.handles.SASetSpecgramOverlapEd,'String',num2str(SASet.Specgram.spect_overlap));
	set(SASet.handles.SASetSpecgramFloorEd,'String',num2str(SASet.Specgram.spect_floor));
	set(SASet.handles.SASetSpecgramCeilEd,'String',num2str(SASet.Specgram.spect_ceil));	
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
	evparamcb('specgram_cmap','set',fig);
	evparamcb('specgram_gamma_style','set',fig);
	set(SASet.handles.SASetSpecgramRangeSld,'Min',0.0,'Max',4.0,'Value',SASet.Specgram.spect_range);
      
      elseif strcmp(varargin{1},'set')
	set(SASet.handles.SASetSpecgramWinDurEd,'String',num2str(SASet.Specgram.spect_win_dur));
	set(SASet.handles.SASetSpecgramOverlapEd,'String',num2str(SASet.Specgram.spect_overlap));
	set(SASet.handles.SASetSpecgramFloorEd,'String',num2str(SASet.Specgram.spect_floor));
	set(SASet.handles.SASetSpecgramCeilEd,'String',num2str(SASet.Specgram.spect_ceil));
	evparamcb('specgram_cmap','set',fig);
	evparamcb('specgram_gamma_style','set',fig);
	set(SASet.handles.SASetSpecgramRangeSld,'Value',SASet.Specgram.spect_range);
      
      elseif strcmp(varargin{1},'get')
	SASet.Specgram.spect_win_dur = str2num(get(SASet.handles.SASetSpecgramWinDurEd,'String'));
	SASet.Specgram.spect_overlap = str2num(get(SASet.handles.SASetSpecgramOverlapEd,'String'));
	SASet.Specgram.spect_floor = str2num(get(SASet.handles.SASetSpecgramFloorEd,'String'));
	SASet.Specgram.spect_ceil = str2num(get(SASet.handles.SASetSpecgramCeilEd,'String'));
	SASet.Specgram.spect_cmap_name = SASet.Specgram.spect_cmap_name_strs{get(SASet.handles.SASetSpecgramCmapPu,'Value')};
	SASet.Specgram.spect_gamma_type = SASet.Specgram.spect_gamma_type_strs{get(SASet.handles.SASetSpecgramGammaStylePu,'Value')};
	SASet.Specgram.spect_range = get(SASet.handles.SASetSpecgramRangeSld,'Value');
	% The 'get' method updates the appdata.
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      end
    end

  case 'specgram_overlap'
    if nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');
    SASet.Specgram.spect_overlap = str2num(get(SASet.handles.SASetSpecgramOverlapEd,'String'));
    if SASet.Specgram.spect_overlap < 0.0
      SASet.Specgram.spect_overlap = 0.0;
    elseif SASet.Specgram.spect_overlap >= 1.0
      SASet.Specgram.spect_overlap = 0.99;
    end
    set(SASet.handles.SASetSpecgramOverlapEd,'String',num2str(SASet.Specgram.spect_overlap));
    setappdata(SASet.handles.SASet,'SASet_Data',SASet);
    
  case 'specgram_floor'
    if nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');
    SASet.Specgram.spect_floor = str2num(get(SASet.handles.SASetSpecgramFloorEd,'String'));
    if SASet.Specgram.spect_floor < -100
      SASet.Specgram.spect_floor = -100;
    elseif SASet.Specgram.spect_floor > 0
      SASet.Specgram.spect_floor = 0;
    end
    set(SASet.handles.SASetSpecgramFloorEd,'String',num2str(SASet.Specgram.spect_floor));
    setappdata(SASet.handles.SASet,'SASet_Data',SASet);
    % Make sure ceiling value is consistent
    evparamcb('specgram_ceil');
    
  case 'specgram_ceil'
    if nargin == 2
      fig = varargin{1};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');
    SASet.Specgram.spect_ceil = str2num(get(SASet.handles.SASetSpecgramCeilEd,'String'));
    if SASet.Specgram.spect_ceil < SASet.Specgram.spect_floor
      SASet.Specgram.spect_ceil = SASet.Specgram.spect_floor;
    elseif SASet.Specgram.spect_ceil > 0
      SASet.Specgram.spect_ceil = 0;
    end
    set(SASet.handles.SASetSpecgramCeilEd,'String',num2str(SASet.Specgram.spect_ceil));
    setappdata(SASet.handles.SASet,'SASet_Data',SASet);
  
  case 'specgram_cmap'
    if nargin == 2
      fig = varargin{1};
    elseif nargin == 3
      myaction = varargin{1}
      fig = varargin{2};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');

    if exist('myaction','var') & strcmp(myaction,'set')
      cmapidx = find(strcmpi(SASet.Specgram.spect_cmap_name_strs,SASet.Specgram.spect_cmap_name));
      if ~isempty(cmapidx)
	set(SASet.handles.SASetSpecgramCmapPu,'Value',cmapidx);
      end
    end
 
    % If colormap is not 'gray' then set the gamma style and
    % disable the gamma popup
    SASet.Specgram.spect_cmap_name = SASet.Specgram.spect_cmap_name_strs{get(SASet.handles.SASetSpecgramCmapPu,'Value')};
    if ~strcmpi(popupstr(SASet.handles.SASetSpecgramCmapPu),'Gray')
      SASet.Specgram.spect_gamma_type = 'brightness';
      setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      evparamcb('specgram_gamma_style','set',fig);
      set(SASet.handles.SASetSpecgramGammaStylePu,'Enable','off');
    else
      set(SASet.handles.SASetSpecgramGammaStylePu,'Enable','on');
    end
      
  case 'specgram_gamma_style'
    if nargin == 2
      fig = varargin{1};
    elseif nargin == 3
      myaction = varargin{1}
      fig = varargin{2};
    else
      fig = gcbf;
    end
    SASet = getappdata(fig,'SASet_Data');

    if exist('myaction','var') & strcmp(myaction,'set')
      gammaidx = find(strcmpi(SASet.Specgram.spect_gamma_type_strs,SASet.Specgram.spect_gamma_type));
      if ~isempty(gammaidx)
	set(SASet.handles.SASetSpecgramGammaStylePu,'Value',gammaidx);
      end
    end
    
    % Depending on the gamma correction type, set fake range on the
    % correction exponent slider.
    SASet.Specgram.spect_gamma_type = SASet.Specgram.spect_gamma_type_strs{get(SASet.handles.SASetSpecgramGammaStylePu,'Value')};
    if strcmpi(SASet.Specgram.spect_gamma_type,'classic')
      set(SASet.handles.SASetSpecgramGammaSldMinTxt,'String',' 0');
      set(SASet.handles.SASetSpecgramGammaSldMaxTxt,'String','4');
      set(SASet.handles.SASetSpecgramExpTxt,'String','Correction exponent');
    elseif strcmpi(SASet.Specgram.spect_gamma_type,'brightness')
      set(SASet.handles.SASetSpecgramGammaSldMinTxt,'String','-1');
      set(SASet.handles.SASetSpecgramGammaSldMaxTxt,'String','1');
      set(SASet.handles.SASetSpecgramExpTxt,'String','Beta value');
    end
  
  % Printing parameters
  case 'print_params'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'init')
	SASet.Print.pr_orient = SASet.params.pr_orient;
	SASet.Print.pr_papertype = SASet.params.pr_papertype;
	SASet.Print.pr_y_paperinches = SASet.params.pr_y_paperinches;
	SASet.Print.pr_x_per_paperinch = SASet.params.pr_x_per_paperinch;
	SASet.Print.pr_left_margin = SASet.params.pr_left_margin;
	SASet.Print.pr_bottom_margin = SASet.params.pr_bottom_margin;
	SASet.Print.pr_orient_strs = {'landscape', 'portrait'};
	SASet.Print.pr_papertype_strs = {'usletter', 'uslegal', 'a4'};	
	set(SASet.handles.SASetPrintYPaperInEd,'String',num2str(SASet.Print.pr_y_paperinches));
	set(SASet.handles.SASetPrintXPerPaperInEd,'String',num2str(SASet.Print.pr_x_per_paperinch));
	set(SASet.handles.SASetPrintLeftMarginEd,'String',num2str(SASet.Print.pr_left_margin));
	set(SASet.handles.SASetPrintBottomMarginEd,'String',num2str(SASet.Print.pr_bottom_margin));

	ptidx = find(strcmpi(SASet.Print.pr_papertype_strs,SASet.Print.pr_papertype));
	if ~isempty(ptidx)
	  set(SASet.handles.SASetPrintPaperTypePu,'Value',ptidx);
	end
	oridx = find(strcmpi(SASet.Print.pr_orient_strs,SASet.Print.pr_orient));
	if ~isempty(oridx)
	  set(SASet.handles.SASetPrintOrientPu,'Value',oridx);
	end
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      
      elseif strcmp(varargin{1},'set')
	set(SASet.handles.SASetPrintYPaperInEd,'String',num2str(SASet.Print.pr_y_paperinches));
	set(SASet.handles.SASetPrintXPerPaperInEd,'String',num2str(SASet.Print.pr_x_per_paperinch));
	set(SASet.handles.SASetPrintLeftMarginEd,'String',num2str(SASet.Print.pr_left_margin));
	set(SASet.handles.SASetPrintBottomMarginEd,'String',num2str(SASet.Print.pr_bottom_margin));

	ptidx = find(strcmpi(SASet.Print.pr_papertype_strs,SASet.Print.pr_papertype));
	if ~isempty(ptidx)
	  set(SASet.handles.SASetPrintPaperTypePu,'Value',ptidx);
	end
	oridx = find(strcmpi(SASet.Print.pr_orient_strs,SASet.Print.pr_orient));
	if ~isempty(oridx)
	  set(SASet.handles.SASetPrintOrientPu,'Value',oridx);
	end
      
      elseif strcmp(varargin{1},'get')
	SASet.Print.pr_y_paperinches = str2num(get(SASet.handles.SASetPrintYPaperInEd,'String'));
	SASet.Print.pr_x_per_paperinch = str2num(get(SASet.handles.SASetPrintXPerPaperInEd,'String'));
	SASet.Print.pr_left_margin = str2num(get(SASet.handles.SASetPrintLeftMarginEd,'String'));
	SASet.Print.pr_bottom_margin = str2num(get(SASet.handles.SASetPrintBottomMarginEd,'String'));
	SASet.Print.pr_papertype = SASet.Print.pr_papertype_strs{get(SASet.handles.SASetPrintPaperTypePu,'Value')};
	SASet.Print.pr_orient = SASet.Print.pr_orient_strs{get(SASet.handles.SASetPrintOrientPu,'Value')};
	% The 'get' method updates the appdata.
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      end
    end

  % Miscellaneous parameters
  case 'misc_params'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'init')
	SASet.Misc.d_Fs = SASet.params.d_Fs;
	SASet.Misc.xrange_start = SASet.params.xrange_start;
	SASet.Misc.do_zc = SASet.params.do_zc;
	SASet.Misc.zcwin = SASet.params.zcwin;	
	set(SASet.handles.SASetMiscDFsEd,'String',num2str(SASet.Misc.d_Fs));
	set(SASet.handles.SASetMiscXrangeEd,'String',num2str(SASet.Misc.xrange_start),'Enable', 'off', 'BackgroundColor', [0.8,0.8,0.8]);
	set(SASet.handles.SASetMiscDoZCCkb,'Value',SASet.Misc.do_zc);
	set(SASet.handles.SASetMiscZCWinEd,'String',num2str(SASet.Misc.zcwin));
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
	evparamcb('do_zc', 'set', fig);
	
      elseif strcmp(varargin{1},'set')
	set(SASet.handles.SASetMiscDFsEd,'String',num2str(SASet.Misc.d_Fs));
	set(SASet.handles.SASetMiscXrangeEd,'String',num2str(SASet.Misc.xrange_start));
	set(SASet.handles.SASetMiscDoZCCkb,'Value',SASet.Misc.do_zc);
	set(SASet.handles.SASetMiscZCWinEd,'String',num2str(SASet.Misc.zcwin));
	evparamcb('do_zc', 'set', fig);
      
      elseif strcmp(varargin{1},'get')
	SASet.Misc.d_Fs = str2num(get(SASet.handles.SASetMiscDFsEd,'String'));
	SASet.Misc.xrange_start = str2num(get(SASet.handles.SASetMiscXrangeEd,'String'));
	SASet.Misc.do_zc = get(SASet.handles.SASetMiscDoZCCkb,'Value');
	SASet.Misc.zcwin = str2num(get(SASet.handles.SASetMiscZCWinEd,'String'));
	% The 'get' method updates the appdata.
	setappdata(SASet.handles.SASet,'SASet_Data',SASet);
      end
    end

  case 'do_zc'
    if nargin == 3
      fig = varargin{2};
      SASet = getappdata(fig,'SASet_Data');
      if strcmp(varargin{1},'set')
	do_zc = SASet.Misc.do_zc;
      end
    else
      SASet = getappdata(gcbf,'SASet_Data');
      do_zc = get(SASet.handles.SASetMiscDoZCCkb,'Value');
    end
  
    if do_zc == 0
      set(SASet.handles.SASetMiscZCWinEd,'Enable','off','BackgroundColor',[0.8,0.8,0.8]);
    else
      set(SASet.handles.SASetMiscZCWinEd,'Enable','on','BackgroundColor',[1.0,1.0,1.0]);
    end  
    
  case 'create_tabs'
    % Idea here is that given the sheet positions and figure,
    % determine and save layout of tab buttons and selector
    
    fig = varargin{1};
    defaultTabNum = varargin{2};
    callback = 'evparamcb';

    SASet = getappdata(fig,'SASet_Data');
    % [fig, sheetpos] = tabdlg('create', tabstrings, tabdims, ...
    %	'tabcbk', [600, 425], [25,10,25,100],1) 
    SASet.tabdims = i_DetermineTabDims(SASet.tabstrings);    
    hsheet = findobj(fig,'Tag','SASetSheet');
    i_CreateTabs(SASet.tabstrings, SASet.tabdims, callback, defaultTabNum, ...
	fig, hsheet);

    setappdata(fig,'SASet_Data',SASet);      
  
  case 'reset_tabdlg'
    fig = varargin{1};
    defaultTab = varargin{2};
    defaultTabNum = varargin{3};
    % DialogUserData = i_GetDialogData(fig);
    % DialogUserData.callback = 'evparamcb';
    % i_SetDialogData(fig, DialogUserData);
    tabdlg('tabpress',fig, defaultTab, defaultTabNum);
    
  case 'tabcallbk'
    pressedTab = varargin{1};
    pressedTabNum = varargin{2};
    previousTab = varargin{3};
    previousTabNum = varargin{4};
    fig = varargin{5};    
    SASet = getappdata(fig,'SASet_Data');

    ntabs = length(SASet.tabfields);
    if previousTabNum > ntabs | pressedTabNum > ntabs
      warn('Tab number out of range.')
      return;
    end
    sprevhandles = ['SASet.', SASet.tabfields{previousTabNum}, '.handles']; 
    scurrhandles = ['SASet.', SASet.tabfields{pressedTabNum}, '.handles'];

    % Turn off prev page's uicontrols and turn on the new one's.
    if isfield(SASet, SASet.tabfields{previousTabNum})
      eval(['prevhandles = ', sprevhandles,';']); 
      if ~isempty(prevhandles)
	set(prevhandles,'Visible','off');
      end
    end
    if isfield(SASet, SASet.tabfields{pressedTabNum})
      eval(['currhandles = ', scurrhandles,';']); 
      if ~isempty(currhandles)
	set(currhandles,'Visible','on');
      end
    end

  otherwise
   error('Unknown file callback action.')
   
end

%
% Helper stuff from tabdlg.m:
%

%******************************************************************************
% Function - Get the user data for the tabbed dialog.                       ***
%******************************************************************************
function data = i_GetDialogData(dialog),

oldHiddenHandleStatus = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');

dataContainer = findobj(dialog,...
  'Type',       'uicontrol', ...
  'Style',      'text', ...
  'Tag',        'TMWDlgDat@#' ...
);

data = get(dataContainer, 'UserData');

set(0, 'ShowHiddenHandles', oldHiddenHandleStatus);


%******************************************************************************
% Function - Set the user data for the tabbed dialog.                       ***
%******************************************************************************
function i_SetDialogData(dialog, data),

oldHiddenHandleStatus = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');

dataContainer = findobj(dialog,...
  'Type',       'uicontrol', ...
  'Style',      'text', ...
  'Tag',        'TMWDlgDat@#' ...
);

if isempty(dataContainer),
  dataContainer = uicontrol(...
    'Parent',           dialog, ...
    'Style',            'text', ...
    'Visible',          'off', ...
    'Tag',              'TMWDlgDat@#' ...
  );
end

set(dataContainer, 'UserData', data);

set(0, 'ShowHiddenHandles', oldHiddenHandleStatus);

%******************************************************************************
% Function - Press the specified tab.                                       ***
%******************************************************************************
function DialogUserData = i_PressTab(hfig, DialogUserData, pressedTabNum),

tabunits_orig = get(DialogUserData.tabs(pressedTabNum), 'Units');
selunits_orig = get(DialogUserData.selector, 'Units');
set(DialogUserData.tabs(pressedTabNum), 'Units', 'pixels');
set(DialogUserData.selector, 'Units', 'pixels');

posPressedTab = get(DialogUserData.tabs(pressedTabNum), 'Position');

posPressedTab(1) = posPressedTab(1) - DialogUserData.selectionHoffset;
posPressedTab(3) = posPressedTab(3) + DialogUserData.selectionHoffset;
posPressedTab(4) = posPressedTab(4) + DialogUserData.selectionVoffset;

set(DialogUserData.tabs(pressedTabNum), 'Position', posPressedTab);

set(DialogUserData.selector, ...
  'Position',           DialogUserData.selectorPos(pressedTabNum,:) ...
);

DialogUserData.activeTabNum = pressedTabNum;

set(DialogUserData.tabs(pressedTabNum), 'Units', tabunits_orig);
set(DialogUserData.selector, 'Units', selunits_orig);


%******************************************************************************
% Function - Unpress the specified tab.                                     ***
%                                                                           ***
% Reduces the size of the specified tab.                                    ***
%                                                                           ***
% NOTE: This function does not move the selector or update the              ***
%   activeTabNum field.  It is assumed that a call to i_PressTab will       ***
%   soon occur and take care of these tasks.                                ***
%******************************************************************************
function i_UnPressTab(hTab, nTab, DialogUserData),

tabunits_orig = get(DialogUserData.tabs(nTab), 'Units');
set(DialogUserData.tabs(nTab), 'Units', 'pixels');

posTab = get(DialogUserData.tabs(nTab), 'Position');

posTab(1) = posTab(1) + DialogUserData.selectionHoffset;
posTab(3) = posTab(3) - DialogUserData.selectionHoffset;
posTab(4) = posTab(4) - DialogUserData.selectionVoffset;

set(DialogUserData.tabs(nTab), 'Position', posTab);

set(DialogUserData.tabs(nTab), 'Units', tabunits_orig);

%******************************************************************************
% Function - Process tab press action.                                      ***
%******************************************************************************
function [DialogUserData, bModified] = ...
  i_ProcessTabPress(hfig, DialogUserData, string, pressedTabNum),

%==============================================================================
% Initialize.
%==============================================================================
bModified = 0;

tabs         = DialogUserData.tabs;
activeTabNum = DialogUserData.activeTabNum;

if pressedTabNum == activeTabNum,
  return;
end

i_UnPressTab(tabs(activeTabNum), activeTabNum, DialogUserData);
DialogUserData = i_PressTab(hfig, DialogUserData, pressedTabNum);
bModified = 1;

%******************************************************************************
% Function - Determine the widths of the tabs based on the strings.         ***
%******************************************************************************
function tabdims = i_DetermineTabDims(strings, font),

%==============================================================================
% Argument checks.
%==============================================================================
if nargin == 1,
  fontsize = get(0, 'DefaultUicontrolFontSize');
  fontname = get(0, 'DefaultUicontrolFontName');
else
  fontsize = font{2};
  fontname = font{1};
end

%==============================================================================
% Create figure and sample text control.
%==============================================================================
hfig  = figure('Visible', 'off');
hText = uicontrol('Style', 'text', 'FontWeight', 'bold');

%==============================================================================
% Get widths.
%==============================================================================
tabdims{1} = zeros(length(strings), 1);
for i=1:length(strings),
  set(hText, 'String', strings{i});
  ext = get(hText, 'Extent');
  tabdims{1}(i) = ext(3) + 4;
end  

tabdims{2} = ext(4) + 2;

%==============================================================================
% Delete objects.
%==============================================================================
delete(hfig);





%******************************************************************************
% Function - Create the tabs for an existing dialog box.
%******************************************************************************
function i_CreateTabs(strings, tabDims, callback, default_page, ...
    hfig, hsheet, font)

%==============================================================================
% Argument checks.
%==============================================================================
if nargin >= 7 & ~isempty(font),
  fontsize = font{2};
  fontname = font{1};
else
  fontsize = get(0, 'FactoryUicontrolFontSize');
  fontname = get(0, 'FactoryUicontrolFontName');
end

origFigUnits = get(hfig, 'Units');
origDefaultUicontrolEnable = get(0, 'DefaultUicontrolEnable');
origDefaultUicontrolUnits = get(0, 'DefaultUicontrolUnits');

set(hfig, ...
    'Units',                              'pixels', ...
    'DefaultUicontrolUnits',              'pixels', ...
    'DefaultUicontrolEnable',             'inactive' ... 
    );

%==============================================================================
% Calculate geometry constants.
%==============================================================================
stringHeight  = tabDims{2};
tabHeight     = tabDims{2} + 5;
tabWidths     = [0; tabDims{1}(:)];
numTabs       = length(tabWidths) - 1;


switch(computer),

  case 'PCWIN',
    
    leftBevelOffset         = 0;
    rightBevelOffset        = 2;   
    topBevelOffset          = 1;
%    selectorHeight          = 4;
    selectorHeight          = 2;
    selectorLeftFudgeFactor = 0;
    deltaTabs               = 3;
    selectionHoffset        = 1;
    sheetEnableState        = 'off';

  case 'MAC2',
   
    leftBevelOffset         = 2;
    rightBevelOffset        = 2;
    topBevelOffset          = 2;
    selectorHeight          = 3;
    selectorLeftFudgeFactor = -1;
    deltaTabs               = 1;
    selectionHoffset        = deltaTabs;
    sheetEnableState        = 'inactive';
 
  otherwise,

    leftBevelOffset         = 3;
    rightBevelOffset        = 3;
    topBevelOffset          = 2;
%    topBevelOffset          = 2;
%    selectorHeight          = 4;
    selectorHeight          = 2;
    selectorLeftFudgeFactor = 0;
    deltaTabs               = 0;
    selectionHoffset        = deltaTabs;
    sheetEnableState        = 'off';

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to give the selected tab a 3-D look, it is made slightly 
%  taller & wider than its unselected size.  selectionVoffset is the
%  number of pixels by which the selected tabs height is increased.
%  Likewise for selectionHoffset.
% NOTE: The 1st tab only lines up w/ the left side of the sheet when
%       it is the selected tab!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectionVoffset = 2;

% Sheet dimensions
origSheetUnits = get(hsheet, 'Units');
set(hsheet, 'Units', 'pixels');
sheetPos = get(hsheet, 'Position');

sheetPos(4) = sheetPos(4) - topBevelOffset;
sheetDims(1) = sheetPos(3);
sheetDims(2) = sheetPos(4);

%==============================================================================
% Create the tabs & store selector positions.
%==============================================================================
posTab(4) = tabHeight;
posTab(2) = sheetPos(2) + sheetPos(4) - 1;

tabs(numTabs) = 0;
selectorPos   = zeros(numTabs, 4);

for i = 1:numTabs,

  butDownFcn = ['tabdlg(''tabpress'', gcbf, ''' strings{i} ''', ' sprintf('%d',i) ')'];

  leftEdge =...
    sheetPos(1)            + ...
    selectionHoffset       + ...
    sum(tabWidths(1:i))    + ...
    ( (i-1) * deltaTabs );

  posTab(1) = leftEdge;
  posTab(3) = tabWidths(i+1);

  tabs(i) = uicontrol( ...
    'Parent',           hfig, ...
    'String',           strings{i}, ...
    'Position',         posTab, ...
    'HorizontalAlign',  'center', ...
    'ButtonDownFcn',    butDownFcn ...
  );

  selectorPos(i, :) = [ ...
    leftEdge - selectionHoffset + leftBevelOffset + selectorLeftFudgeFactor ...
    posTab(2) ...
    posTab(3) + selectionHoffset - rightBevelOffset ...
    selectorHeight ...
  ];

end

% Put the tab buttons behind everything (matlab6.0 has a bug in uistack!)
if str2num(version('-release')) >= 12.1  
  uistack(tabs,'bottom');
end

%==============================================================================
% Set sheet defaults.
%==============================================================================
set(hsheet, 'Enable', sheetEnableState);

%==============================================================================
% Create the selector.
%==============================================================================
selector = uicontrol( ...
  'Parent',             hfig, ...
  'Style',              'text' ...
);

%==============================================================================
% Save pertinent info in tabbed dialog data container.
%==============================================================================
DialogUserData.tabs             = tabs;
DialogUserData.selector         = selector;
DialogUserData.selectorPos      = selectorPos;
DialogUserData.selectionHoffset = selectionHoffset;
DialogUserData.selectionVoffset = selectionVoffset;
DialogUserData.leftBevelOffset  = leftBevelOffset;
DialogUserData.rightBevelOffset = rightBevelOffset;
DialogUserData.deltaTabs        = deltaTabs;
DialogUserData.activeTabNum     = -1;
DialogUserData.callback         = callback;
DialogUserData.strings          = strings;

%==============================================================================
% Select the default tab.
%==============================================================================
DialogUserData = i_PressTab(hfig, DialogUserData, default_page);

%==============================================================================
% Store the user data.
%==============================================================================
i_SetDialogData(hfig, DialogUserData);

%==============================================================================
% Restore defaults.
%==============================================================================
set(hfig, 'DefaultUicontrolEnable', origDefaultUicontrolEnable);
set(hfig, 'DefaultUicontrolUnits', origDefaultUicontrolUnits);

% Restore state
set(hfig, 'Units', origFigUnits);

% Set everything to characters units after set up.
set(tabs, 'Units', 'characters');
set(hsheet, 'Units', origSheetUnits);
set(selector, 'Units', 'characters');


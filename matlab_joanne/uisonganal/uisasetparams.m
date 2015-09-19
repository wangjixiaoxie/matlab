function varargout = uisasetparams(varargin)
% UISASETPARAMS Application M-file for uisasetparams.fig
%    FIG = UISASETPARAMS launch uisasetparams GUI.
%    UISASETPARAMS('callback_name', ...) invoke the named callback.

% This function was created by saveas(...,'mfig'), or print -dmfile.
% It loads an HG object saved in binary format in the FIG-file of the
% same name.  NOTE: if you want to see the old-style M code
% representation of a saved object, previously created by print -dmfile,
% you can obtain this by using saveas(...,'mmat'). But be advised that the
% M-file/MAT-file format does not preserve some important information due
% to limitations in that format, including ApplicationData stored using
% setappdata.  Also, references to handles stored in UserData or Application-
% Data will no longer be valid if saved in the M-file/MAT-file format.
% ApplicationData and stored handles (excluding integer figure handles)
% are both correctly saved in FIG-files.
%

if nargin == 0  % LAUNCH GUI

  %load the saved object
  [path, name] = fileparts(which(mfilename));
  figname = fullfile(path, [name '.fig']);
  if (exist(figname,'file'))
    fig = openfig(figname,'reuse') 
  else 
    fig = openfig([name '.fig'],'reuse')
  end
  % Make sure the dialog is on screen
  movegui(fig,'onscreen')
  
  % Use system color scheme for figure:
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

  % Generate a structure of handles to pass to callbacks, and store it. 
  SASet.handles = guihandles(fig);  
  setappdata(fig, 'SASet_Data', SASet);

  if nargout > 0
    varargout{1} = fig;
  end

  % Wait for callbacks to run and window to be dismissed:
  % uiwait(fig);

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

  try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
  catch
    disp(lasterr);
  end
  
end

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox5_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox6_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox8_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramWinDurEd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramOverlapEd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramFloorEd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramCeilingEd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramCmapPu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramGammaStylePu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetSpecgramRangeSld_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetMiscDFsEd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = SASetMiscZCWinEd_Callback(h, eventdata, handles, varargin)


function varargout = uifig(varargin)
% UIFIG Application M-file for uifig.fig
%    FIG = UIFIG launch uifig GUI.
%    UIFIG('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 07-May-2001 22:51:13

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'new');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

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
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton3.
disp('pushbutton3 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton4.
disp('pushbutton4 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton5.
disp('pushbutton5 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = pushbutton6_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton6.
disp('pushbutton6 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = pushbutton7_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton7.
disp('pushbutton7 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton1.
disp('radiobutton1 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton2.
disp('radiobutton2 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton4.
disp('radiobutton4 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton5.
disp('radiobutton5 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton6_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton6.
disp('radiobutton6 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton7_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton7.
disp('radiobutton7 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton8_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton8.
disp('radiobutton8 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton9_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton9.
disp('radiobutton9 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton10_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton10.
disp('radiobutton10 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton11_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton11.
disp('radiobutton11 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton12_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton12.
disp('radiobutton12 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton13_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton13.
disp('radiobutton13 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = radiobutton14_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton14.
disp('radiobutton14 Callback not implemented yet.')


% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu1.
disp('popupmenu1 Callback not implemented yet.')
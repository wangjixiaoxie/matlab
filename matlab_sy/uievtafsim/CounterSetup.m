function varargout = CounterSetup(varargin)
% COUNTERSETUP M-file for CounterSetup.fig
%      COUNTERSETUP, by itself, creates a new COUNTERSETUP or raises the existing
%      singleton*.
%
%      H = COUNTERSETUP returns the handle to a new COUNTERSETUP or the handle to
%      the existing singleton*.
%
%      COUNTERSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COUNTERSETUP.M with the given input arguments.
%
%      COUNTERSETUP('Property','Value',...) creates a new COUNTERSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CounterSetup_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CounterSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CounterSetup

% Last Modified by GUIDE v2.5 18-Oct-2006 18:49:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CounterSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @CounterSetup_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CounterSetup is made visible.
function CounterSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CounterSetup (see VARARGIN)

CntStruct=varargin{1};
handles.OrigCntRng=CntStruct; %counter range strucutre
handles.CntRng=CntStruct;
handles.CurInd=1;
handles.output=CntStruct;
set(handles.WhichText,'String',['1/',num2str(length(CntStruct))]);
guidata(hObject, handles);
SetVals(handles);

if (handles.CurInd==length(handles.CntRng))
    set(handles.DOAND,'Enable','off');
end

uiwait(handles.CounterSetup);

% UIWAIT makes CounterSetup wait for user response (see UIRESUME)
%uiwait(handles.CounterSetup);
return;

function SetVals(handles);
% sets the current values on the screen

ind=handles.CurInd;
CntStruct=handles.CntRng;
set(handles.MinCntBox,'String',num2str(CntStruct(ind).MinCnt));
set(handles.MaxCntBox,'String',num2str(CntStruct(ind).MaxCnt));
set(handles.ThresholdBox,'String',num2str(CntStruct(ind).Thresh));
if (CntStruct(ind).Mode==1)
    set(handles.ModeBox,'Value',get(handles.ModeBox,'Min'));
else
    set(handles.ModeBox,'Value',get(handles.ModeBox,'Max'));
end

if (CntStruct(ind).DoNot==1)
    set(handles.NotBox,'Value',get(handles.NotBox,'Max'));
else
    set(handles.NotBox,'Value',get(handles.NotBox,'Min'));
end

if (CntStruct(ind).DoAND==1)
    set(handles.DOAND,'Value',get(handles.DOAND,'Max'));
else
    set(handles.DOAND,'Value',get(handles.DOAND,'Min'));
end
set(handles.WhichText,'String',[num2str(ind),'/',num2str(length(CntStruct))]);
return;


% --- Outputs from this function are returned to the command line.
function varargout = CounterSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.CounterSetup);
return;

% --- Executes during object creation, after setting all properties.
function WhichText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhichText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
return;

% --- Executes on button press in PrevBtn.
function PrevBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PrevBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.CurInd>1)
    handles.CurInd=handles.CurInd - 1;
    guidata(hObject,handles);
    SetVals(handles);

    if (handles.CurInd==length(handles.CntRng))
        set(handles.DOAND,'Enable','off');
    else
        set(handles.DOAND,'Enable','on');
    end
end
return;

% --- Executes on button press in NextBtn.
function NextBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NextBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.CurInd<length(handles.CntRng))
    handles.CurInd=handles.CurInd + 1;
    guidata(hObject,handles);
    SetVals(handles);

    if (handles.CurInd==length(handles.CntRng))
        set(handles.DOAND,'Enable','off');
    else
        set(handles.DOAND,'Enable','on');
    end
end
return;

% --- Executes on button press in NotBox.
function NotBox_Callback(hObject, eventdata, handles)
% hObject    handle to NotBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NotBox

ind=handles.CurInd;
if (get(handles.NotBox,'Value')==get(handles.NotBox,'Max'))
    handles.CntRng(ind).DoNot=1;
else
    handles.CntRng(ind).DoNot=0;
end
guidata(hObject,handles);
return;

% --- Executes on button press in ModeBox.
function ModeBox_Callback(hObject, eventdata, handles)
% hObject    handle to ModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of ModeBox

ind=handles.CurInd;
if (get(handles.ModeBox,'Value')==get(handles.ModeBox,'Max'))
    handles.CntRng(ind).Mode=0;
else
    handles.CntRng(ind).Mode=1;
end
guidata(hObject,handles);

return

% --- Executes on button press in OkayBtn.
function OkayBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OkayBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp=handles.CntRng;
handles.output=tmp;
guidata(hObject,handles);
uiresume(handles.CounterSetup);
return;

% --- Executes on button press in CancelBtn.
function CancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp=handles.OrigCntRng;
handles.output=tmp;
guidata(hObject,handles);
uiresume(handles.CounterSetup);
return;


% --- Executes during object creation, after setting all properties.
function MaxCntBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxCntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

% --- Executes during object creation, after setting all properties.
function ThresholdBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

% --- Executes during object creation, after setting all properties.
function MinCntBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinCntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

function MinCntBox_Callback(hObject, eventdata, handles)
% hObject    handle to MinCntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinCntBox as text
%        str2double(get(hObject,'String')) returns contents of MinCntBox as a double

ind=handles.CurInd;
handles.CntRng(ind).MinCnt=str2num(get(handles.MinCntBox,'String'));
guidata(hObject,handles);
return;

function MaxCntBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxCntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxCntBox as text
%        str2double(get(hObject,'String')) returns contents of MaxCntBox as a double

ind=handles.CurInd;
handles.CntRng(ind).MaxCnt=str2num(get(handles.MaxCntBox,'String'));
guidata(hObject,handles);
return;

function ThresholdBox_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdBox as text
%        str2double(get(hObject,'String')) returns contents of ThresholdBox as a double

ind=handles.CurInd;
handles.CntRng(ind).Thresh=str2num(get(handles.ThresholdBox,'String'));
guidata(hObject,handles);

return;


% --- Executes on button press in DOAND.
function DOAND_Callback(hObject, eventdata, handles)
% hObject    handle to DOAND (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DOAND

ind=handles.CurInd;
if (get(handles.DOAND,'Value')==get(handles.DOAND,'Max'))
    handles.CntRng(ind).DoAND=1;
else
    handles.CntRng(ind).DoAND=0;
end
guidata(hObject,handles);
return;

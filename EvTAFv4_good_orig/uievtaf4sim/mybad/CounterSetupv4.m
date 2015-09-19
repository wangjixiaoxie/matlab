function varargout = CounterSetupv4(varargin)
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

% Last Modified by GUIDE v2.5 16-Aug-2009 21:43:58

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

ND=varargin{1};
handles.OrigND=ND; %counter range strucutre
handles.ND=ND;
handles.CurInd=1;
handles.output=ND;
set(handles.WhichText,'String',['1/',num2str(length(ND.CntRng))]);
guidata(hObject, handles);
SetVals(handles);

uiwait(handles.CounterSetup);

% UIWAIT makes CounterSetup wait for user response (see UIRESUME)
%uiwait(handles.CounterSetup);
return;

function SetVals(handles);
% sets the current values on the screen

ind=handles.CurInd;
CntStruct=handles.ND.CntRng;
set(handles.MinCntBox,'String',num2str(CntStruct(ind).Min));
set(handles.MaxCntBox,'String',num2str(CntStruct(ind).Max));
set(handles.ThresholdBox,'String',num2str(CntStruct(ind).TH));
if (CntStruct(ind).Mode==1)
    set(handles.ModeBox,'Value',get(handles.ModeBox,'Min'));
else
    set(handles.ModeBox,'Value',get(handles.ModeBox,'Max'));
end

if (CntStruct(ind).Not==1)
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

set(hanldes.CntLogBox,'String',ND.CntLog);
set(handles.FreqRngModeBox,'Value',ND.FreqRng.FreqMode+1);
set(handles.MinFreqBox,'String',num2str(ND.FreqRng.MinFreq));
set(handles.MaxFreqBox,'String',num2str(ND.FreqRng.MaxFreq));
set(handles.FerqTHMinBox,'String',num2str(ND.FreqRng.FreqTHMin));
set(handles.FreqTHMaxBox,'String',num2str(ND.FreqRng.FreqTHMax));
set(handles.NBinBox,'String',num2str(ND.FreqRng.NBins));
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
    handles.ND.CntRng(ind).Not=1;
else
    handles.ND.CntRng(ind).Not=0;
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
    handles.ND.CntRng(ind).Mode=0;
else
    handles.ND.CntRng(ind).Mode=1;
end
guidata(hObject,handles);

return

% --- Executes on button press in OkayBtn.
function OkayBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OkayBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp=handles.ND;
handles.output=tmp;
guidata(hObject,handles);
uiresume(handles.CounterSetup);
return;

% --- Executes on button press in CancelBtn.
function CancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to CancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp=handles.OrigND;
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
handles.ND.CntRng(ind).Min=str2num(get(handles.MinCntBox,'String'));
guidata(hObject,handles);
return;

function MaxCntBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxCntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxCntBox as text
%        str2double(get(hObject,'String')) returns contents of MaxCntBox as a double

ind=handles.CurInd;
handles.ND.CntRng(ind).Max=str2num(get(handles.MaxCntBox,'String'));
guidata(hObject,handles);
return;

function ThresholdBox_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdBox as text
%        str2double(get(hObject,'String')) returns contents of ThresholdBox as a double

ind=handles.CurInd;
handles.ND.FreqRng.Thresh=str2num(get(handles.ThresholdBox,'String'));
guidata(hObject,handles);

return;


function BTAFMinBox_Callback(hObject, eventdata, handles)
% hObject    handle to BTAFMinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BTAFMinBox as text
%        str2double(get(hObject,'String')) returns contents of BTAFMinBox as a double
ind=handles.CurInd;
handles.ND.CntRng(ind).BTAFmin=str2num(get(handles.BTAFMinBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function BTAFMinBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BTAFMinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FreqRngModeBox.
function FreqRngModeBox_Callback(hObject, eventdata, handles)
% hObject    handle to FreqRngModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FreqRngModeBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FreqRngModeBox
ind=handles.CurInd;
handles.ND.FreqRng.FreqMode=get(handles.FreqRngModeBox,'Value')-1;
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function FreqRngModeBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqRngModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CntLogBox_Callback(hObject, eventdata, handles)
% hObject    handle to CntLogBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CntLogBox as text
%        str2double(get(hObject,'String')) returns contents of CntLogBox as a double
ind=handles.CurInd;
handles.ND.CntLog=get(handles.CntLogBox,'String');
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function CntLogBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CntLogBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinFreqBox_Callback(hObject, eventdata, handles)
% hObject    handle to MinFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinFreqBox as text
%        str2double(get(hObject,'String')) returns contents of MinFreqBox as a double
ind=handles.CurInd;
handles.ND.FreqRng.MinFreq=str2num(get(handles.MinFreqBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function MinFreqBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxFreqBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxFreqBox as text
%        str2double(get(hObject,'String')) returns contents of MaxFreqBox as a double
ind=handles.CurInd;
handles.ND.FreqRng.MaxFreq=str2num(get(handles.MaxFreqBox,'String'));
guidata(hObject,handles);
return;


% --- Executes during object creation, after setting all properties.
function MaxFreqBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FreqTHMinBox_Callback(hObject, eventdata, handles)
% hObject    handle to FreqTHMinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqTHMinBox as text
%        str2double(get(hObject,'String')) returns contents of FreqTHMinBox as a double
ind=handles.CurInd;
handles.ND.FreqRng.FreqTHMin=str2num(get(handles.FreqTHMinBox,'String'));
guidata(hObject,handles);
return;


% --- Executes during object creation, after setting all properties.
function FreqTHMinBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqTHMinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FreqTHMaxBox_Callback(hObject, eventdata, handles)
% hObject    handle to FreqTHMaxBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqTHMaxBox as text
%        str2double(get(hObject,'String')) returns contents of FreqTHMaxBox as a double
ind=handles.CurInd;
handles.ND.FreqRng.FreqTHMax=str2num(get(handles.FreqTHMaxBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function FreqTHMaxBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqTHMaxBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NBinBox_Callback(hObject, eventdata, handles)
% hObject    handle to NBinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NBinBox as text
%        str2double(get(hObject,'String')) returns contents of NBinBox as a double
ind=handles.CurInd;
handles.ND.FreqRng.NBins=str2num(get(handles.NBinBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function NBinBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NBinBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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

% Last Modified by GUIDE v2.5 20-Aug-2009 14:30:17

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
NDIndex=varargin{2};
tmpstr=[];
for ijk=1:length(ND)
    tmpstr{ijk}=num2str(ijk);
end
set(handles.NDIndexBox,'String',tmpstr);
set(handles.NDIndexBox,'Value',NDIndex);

%ND=ND(NDIndex);

CntStruct=ND(NDIndex).CntRng;
handles.OrigND=ND; %counter range strucutre
handles.OrigNDIndex=NDIndex; %counter range strucutre
handles.ND=ND;
handles.NDIndex=NDIndex;
handles.CurInd=1; % counter range index
handles.output=ND;

NTempl=size(handles.ND(handles.NDIndex).Templ,2);
handles.NTempl=NTempl;

set(handles.WhichText,'String',['1/',num2str(NTempl)]);
handles.WASCANC=0;
guidata(hObject, handles);
SetVals(handles);

uiwait(handles.CounterSetup);

% UIWAIT makes CounterSetup wait for user response (see UIRESUME)
%uiwait(handles.CounterSetup);
return;

function SetVals(handles);
% sets the current values on the screen

ind=handles.CurInd;
ND=handles.ND(handles.NDIndex);
CntStruct=ND.CntRng;
set(handles.MinCntBox,'String',num2str(ND.CntRng(ind).Min));
set(handles.MaxCntBox,'String',num2str(ND.CntRng(ind).Max));
set(handles.ThresholdBox,'String',num2str(ND.CntRng(ind).TH));
if (ND.CntRng(ind).Mode==1)
    set(handles.ModeBox,'Value',get(handles.ModeBox,'Min'));
else
    set(handles.ModeBox,'Value',get(handles.ModeBox,'Max'));
end

if (ND.CntRng(ind).Not==1)
    set(handles.NotBox,'Value',get(handles.NotBox,'Max'));
else
    set(handles.NotBox,'Value',get(handles.NotBox,'Min'));
end

set(handles.CntLogBox,'String',ND.CntLog);
set(handles.FreqModeBox,'Value',ND.FreqRng.FreqMode+1);
set(handles.FFMinFreqBox,'String',num2str(ND.FreqRng.MinFreq));
set(handles.FFMaxFreqBox,'String',num2str(ND.FreqRng.MaxFreq));
set(handles.FreqTHMinBox,'String',num2str(ND.FreqRng.FreqTHMin));
set(handles.FreqTHMaxBox,'String',num2str(ND.FreqRng.FreqTHMax));
set(handles.NBinBox,'String',num2str(ND.FreqRng.NBins));

set(handles.AmpModeBox,'Value',ND.AmpRng.AmpMode+1);
set(handles.AmpMinFreqBox,'String',num2str(ND.AmpRng.MinFreq));
set(handles.AmpMaxFreqBox,'String',num2str(ND.AmpRng.MaxFreq));
set(handles.AmpTHBox,'String',num2str(ND.AmpRng.AmpThresh));

set(handles.MinRepeatBox,'String',num2str(ND.RepCntRng.MinRepeat));
set(handles.MaxRepeatBox,'String',num2str(ND.RepCntRng.MaxRepeat));
set(handles.RepRefracBox,'String',num2str(ND.RepCntRng.RepRefrac));
set(handles.RepResetBox,'String',num2str(ND.RepCntRng.RepReset));

set(handles.DelayToContingBox,'String',num2str(handles.ND(handles.NDIndex).DelayToContingen));
set(handles.TrigRefracBox,'String',num2str(handles.ND(handles.NDIndex).TrigRefrac));

set(handles.WhichText,'String',[num2str(ind),'/',num2str(handles.NTempl)]);
return;


% --- Outputs from this function are returned to the command line.
function varargout = CounterSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.NDIndex;
varargout{3} = handles.WASCANC;
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

if (handles.CurInd<handles.NTempl)
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
    handles.ND(handles.NDIndex).CntRng(ind).Not=1;
else
    handles.ND(handles.NDIndex).CntRng(ind).Not=0;
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
    handles.ND(handles.NDIndex).CntRng(ind).Mode=0;
else
    handles.ND(handles.NDIndex).CntRng(ind).Mode=1;
end
guidata(hObject,handles);

return

% --- Executes on button press in OkayBtn.
function OkayBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OkayBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp=handles.ND;
handles.output=tmp;handles.WASCANC=0;
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
handles.NDIndex=handles.OrigNDIndex;
handles.WASCANC=1;
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
handles.ND(handles.NDIndex).CntRng(ind).Min=str2num(get(handles.MinCntBox,'String'));
guidata(hObject,handles);
return;

function MaxCntBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxCntBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxCntBox as text
%        str2double(get(hObject,'String')) returns contents of MaxCntBox as a double

ind=handles.CurInd;
handles.ND(handles.NDIndex).CntRng(ind).Max=str2num(get(handles.MaxCntBox,'String'));
guidata(hObject,handles);
return;

function ThresholdBox_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdBox as text
%        str2double(get(hObject,'String')) returns contents of ThresholdBox as a double

ind=handles.CurInd;
handles.ND(handles.NDIndex).CntRng(ind).TH=str2num(get(handles.ThresholdBox,'String'));
guidata(hObject,handles);
return;

function CntLogBox_Callback(hObject, eventdata, handles)
% hObject    handle to CntLogBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CntLogBox as text
%        str2double(get(hObject,'String')) returns contents of CntLogBox as a double
handles.ND(handles.NDIndex).CntLog=get(handles.CntLogBox,'String');
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


% --- Executes on selection change in FreqModeBox.
function FreqModeBox_Callback(hObject, eventdata, handles)
% hObject    handle to FreqModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FreqModeBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FreqModeBox
handles.ND(handles.NDIndex).FreqRng.FreqMode=get(handles.FreqModeBox,'Value')-1;
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function FreqModeBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FFMinFreqBox_Callback(hObject, eventdata, handles)
% hObject    handle to FFMinFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFMinFreqBox as text
%        str2double(get(hObject,'String')) returns contents of FFMinFreqBox as a double
handles.ND(handles.NDIndex).FreqRng.MinFreq=str2num(get(handles.FFMinFreqBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function FFMinFreqBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFMinFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FFMaxFreqBox_Callback(hObject, eventdata, handles)
% hObject    handle to FFMaxFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FFMaxFreqBox as text
%        str2double(get(hObject,'String')) returns contents of FFMaxFreqBox as a double
handles.ND(handles.NDIndex).FreqRng.MaxFreq=str2num(get(handles.FFMaxFreqBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function FFMaxFreqBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FFMaxFreqBox (see GCBO)
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
handles.ND(handles.NDIndex).FreqRng.FreqTHMin=str2num(get(handles.FreqTHMinBox,'String'));
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
handles.ND(handles.NDIndex).FreqRng.FreqTHMax=str2num(get(handles.FreqTHMaxBox,'String'));
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
handles.ND(handles.NDIndex).FreqRng.NBins=str2num(get(handles.NBinBox,'String'));
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


% --- Executes on selection change in AmpModeBox.
function AmpModeBox_Callback(hObject, eventdata, handles)
% hObject    handle to AmpModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns AmpModeBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AmpModeBox
handles.ND(handles.NDIndex).AmpRng.NBins=str2num(get(handles.AmpModeBox,'Value')-1);
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function AmpModeBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AmpModeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AmpMinFreqBox_Callback(hObject, eventdata, handles)
% hObject    handle to AmpMinFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AmpMinFreqBox as text
%        str2double(get(hObject,'String')) returns contents of AmpMinFreqBox as a double
handles.ND(handles.NDIndex).AmpRng.MinFreq=str2num(get(handles.AmpMinFreqBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function AmpMinFreqBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AmpMinFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AmpMaxFreqBox_Callback(hObject, eventdata, handles)
% hObject    handle to AmpMaxFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AmpMaxFreqBox as text
%        str2double(get(hObject,'String')) returns contents of AmpMaxFreqBox as a double
handles.ND(handles.NDIndex).AmpRng.MaxFreq=str2num(get(handles.AmpMaxFreqBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function AmpMaxFreqBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AmpMaxFreqBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AmpTHBox_Callback(hObject, eventdata, handles)
% hObject    handle to AmpTHBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AmpTHBox as text
%        str2double(get(hObject,'String')) returns contents of AmpTHBox as a double
handles.ND(handles.NDIndex).AmpRng.AmpThresh=str2num(get(handles.AmpTHBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function AmpTHBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AmpTHBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinRepeatBox_Callback(hObject, eventdata, handles)
% hObject    handle to MinRepeatBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinRepeatBox as text
%        str2double(get(hObject,'String')) returns contents of MinRepeatBox as a double
handles.ND(handles.NDIndex).RepCntRng.MinRepeat=str2num(get(handles.MinRepeatBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function MinRepeatBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinRepeatBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxRepeatBox_Callback(hObject, eventdata, handles)
% hObject    handle to MaxRepeatBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxRepeatBox as text
%        str2double(get(hObject,'String')) returns contents of MaxRepeatBox as a double
handles.ND(handles.NDIndex).RepCntRng.MaxRepeat=str2num(get(handles.MaxRepeatBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function MaxRepeatBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxRepeatBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RepRefracBox_Callback(hObject, eventdata, handles)
% hObject    handle to RepRefracBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RepRefracBox as text
%        str2double(get(hObject,'String')) returns contents of RepRefracBox as a double
handles.ND(handles.NDIndex).RepCntRng.RepRefrac=str2num(get(handles.RepRefracBox,'String'));
guidata(hObject,handles);
return;


% --- Executes during object creation, after setting all properties.
function RepRefracBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RepRefracBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RepResetBox_Callback(hObject, eventdata, handles)
% hObject    handle to RepResetBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RepResetBox as text
%        str2double(get(hObject,'String')) returns contents of RepResetBox as a double
handles.ND(handles.NDIndex).RepCntRng.RepReset=str2num(get(handles.RepResetBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function RepResetBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RepResetBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NDIndexBox.
function NDIndexBox_Callback(hObject, eventdata, handles)
% hObject    handle to NDIndexBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns NDIndexBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        NDIndexBox

ND=handles.ND;
NDIndex=get(hObject,'Value');

CntStruct=ND(NDIndex).CntRng;

handles.CurInd=1; % counter range index

NTempl=size(handles.ND(NDIndex).Templ,2);
handles.NTempl=NTempl;

handles.ND=ND;
handles.NDIndex=NDIndex;
guidata(hObject, handles);
SetVals(handles);

% --- Executes during object creation, after setting all properties.
function NDIndexBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NDIndexBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TrigRefracBox_Callback(hObject, eventdata, handles)
% hObject    handle to TrigRefracBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrigRefracBox as text
%        str2double(get(hObject,'String')) returns contents of TrigRefracBox as a double
handles.ND(handles.NDIndex).TrigRefrac=str2num(get(handles.TrigRefracBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function TrigRefracBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrigRefracBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DelayToContingBox_Callback(hObject, eventdata, handles)
% hObject    handle to DelayToContingBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DelayToContingBox as text
%        str2double(get(hObject,'String')) returns contents of DelayToContingBox as a double
handles.ND(handles.NDIndex).DelayToContingen=str2num(get(handles.DelayToContingBox,'String'));
guidata(hObject,handles);
return;

% --- Executes during object creation, after setting all properties.
function DelayToContingBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DelayToContingBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddNewNoteButton.
function AddNewNoteButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddNewNoteButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% add a new template to the end

disp(['Choose a template file for the new note']);
[fn,pth]=uigetfile('*.*','Template file for new note');
if (fn==0)
    return;
elseif (~exist(fn,'file'))
    return;
end
newtempl=load(fullfile(pth,fn));

ND=handles.ND;
ND(length(ND)+1)=ND(handles.NDIndex);
NDIndex=length(ND);
tmpstr=[];
for ijk=1:length(ND)
    tmpstr{ijk}=num2str(ijk);
end
set(handles.NDIndexBox,'String',tmpstr);
set(handles.NDIndexBox,'Value',NDIndex);
handles.ND=ND;
handles.NDIndex=NDIndex;
guidata(hObject,handles);

ND=handles.ND;
NDIndex=length(ND);

handles.CurInd=1; % counter range index
ND(NDIndex).Templ=newtempl;
ND(NDIndex).TemplFile=['CHECKPATH__',fn];
NTempl=size(ND(handles.NDIndex).Templ,2);
handles.NTempl=NTempl;

CntRng=ND(NDIndex).CntRng;
CntLog=ND(NDIndex).CntLog;
if (length(CntRng)>NTempl)
    CntRng=CntRng(1:NTempl);
    % make simple CntLog all or's
    cntlogtemp='a';
    for ijk=2:NTempl
        cntlogtemp=[cntlogtemp,' + ',char(fix('a')+ijk-1)];
    end
elseif (length(CntRng)<NTempl)
    szorig=length(CntRng);
    for ijk=szorig+1:NTempl
        CntRng(ijk)=CntRng(szorig);
    end
    % make simple CntLog all or's
    cntlogtemp='a';
    for ijk=2:NTempl
        cntlogtemp=[cntlogtemp,' + ',char(fix('a')+ijk-1)];
    end
end
for ijk=1:length(CntRng)
    CntRng(ijk).VarName=char(fix('a')+ijk-1);
end
ND(NDIndex).CntRng=CntRng;
ND(NDIndex).CntLog=cntlogtemp;

handles.ND=ND;
handles.NDIndex=NDIndex;
guidata(hObject,handles);
SetVals(handles);
return;

% --- Executes on button press in DelNoteButton.
function DelNoteButton_Callback(hObject, eventdata, handles)
% hObject    handle to DelNoteButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ND=handles.ND;
NDIndex=handles.NDIndex;
if (length(ND)==1)
    warndlg('you only have 1 note template - can''t erase it');
    return;
end
ND(NDIndex)=[];
if (NDIndex>length(ND))
    NDIndex=length(ND);
end
handles.ND=ND;
handles.NDIndex=NDIndex;
guidata(hObject,handles);

tmpstr=[];
for ijk=1:length(ND)
    tmpstr{ijk}=num2str(ijk);
end
set(handles.NDIndexBox,'String',tmpstr);
set(handles.NDIndexBox,'Value',NDIndex);

SetVals(handles);
return

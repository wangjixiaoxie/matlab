function varargout = uievscreen(varargin)
% uievscreen('batch','obs0');

% UIEVSCREEN M-file for uievscreen.fig
%      UIEVSCREEN, by itself, creates a new UIEVSCREEN or raises the existing
%      singleton*.
%
%      H = UIEVSCREEN returns the handle to a new UIEVSCREEN or the handle to
%      the existing singleton*.
%
%      UIEVSCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UIEVSCREEN.M with the given input arguments.
%
%      UIEVSCREEN('Property','Value',...) creates a new UIEVSCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before uievscreen_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to uievscreen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help uievscreen

% Last Modified by GUIDE v2.5 03-Aug-2005 19:11:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uievscreen_OpeningFcn, ...
                   'gui_OutputFcn',  @uievscreen_OutputFcn, ...
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


% --- Executes just before uievscreen is made visible.
function uievscreen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to uievscreen (see VARARGIN)

% Choose default command line output for uievscreen
handles.output = hObject;

bt=varargin{1};
CS=varargin{2};

fid=fopen(bt,'r');
if (fid==-1)
    disp(['could not find batch file :',bt]);
    return;
end

cnt=0;
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    cnt=cnt+1;
    fstruct(cnt).fname=fn;
end
fclose(fid);

if (cnt<1)
    disp(['no files exist']);
    return;
end

handles.FILES=fstruct;
handles.NFILE=1;
handles.CHANSPEC=CS;
guidata(hObject, handles);

uievscreen_plot(hObject,handles);
% UIWAIT makes uievscreen wait for user response (see UIRESUME)
uiwait(handles.uievscreen);

return;


% --- Outputs from this function are returned to the command line.
function varargout = uievscreen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.uievscreen);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uievscreen_plot(hObject,handles);

set(handles.FileNumBox,'String',[num2str(handles.NFILE),'/',num2str(length(handles.FILES))]);

[dat,fs]=evsoundin('',handles.FILES(handles.NFILE).fname,handles.CHANSPEC);
if (strcmp(handles.CHANSPEC(1:min([3,length(handles.CHANSPEC)])),'obs'))
    rd=readrecf(handles.FILES(handles.NFILE).fname);
else
    rd=[];
end
if (get(handles.ShowSpecBtn,'Value')==get(handles.ShowSpecBtn,'Max'))
    [sm,sp,t,f]=evsmooth(dat,fs,str2num(get(handles.SPTHBox,'String')));
    handles.sm=sm;
    handles.sp=sp;
    handles.t=t;
    handles.f=f;
    axes(handles.SpecGramAxes);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
else
    [sm]=evsmooth(dat,fs,str2num(get(handles.SPTHBox,'String')));
    handles.sm=sm;
end


handles.dat=dat;
handles.fs=fs;
axes(handles.SmoothAxes);
if (~isfield(handles,'SMYLIM'))
    semilogy([1:10:length(sm)]/fs,sm(1:10:end),'b-');
    v=axis;
    handles.SMYLIM=v(3:4);
else
    v=axis;
    handles.SMYLIM=v(3:4);
    semilogy([1:10:length(sm)]/fs,sm(1:10:end),'b-');
    ylim(v(3:4));
end
if (isfield(rd,'tt'))
    tt=rd.ttimes;
else
    tt=[];
end
if (length(tt)>0)
    yy=ylim;hold on;
    plot(tt/1e3,yy(1)+0.1*diff(yy)+0*tt,'r^');
    hold off;
end
tmp=title(handles.FILES(handles.NFILE).fname,'Interpreter','none');
zoom xon;
linkaxes([handles.SmoothAxes,handles.SpecGramAxes],'x');
guidata(hObject,handles);
return;

% --- Executes on button press in NextFileBtn.
function NextFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.NFILE<length(handles.FILES))
    handles.NFILE=handles.NFILE+1;
    guidata(hObject,handles);
    uievscreen_plot(hObject,handles);
end
return;

% --- Executes on button press in PrevFileBtn.
function PrevFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.NFILE>1)
    handles.NFILE=handles.NFILE-1;
    guidata(hObject,handles);
    uievscreen_plot(hObject,handles);
end
return;

% --- Executes on button press in SkipToBtn.
function SkipToBtn_Callback(hObject, eventdata, handles)
% hObject    handle to SkipToBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[tmpout,isCancel] = SkipToFile(handles.FILES,handles.NFILE);
if (~isCancel)
    handles.NFILE=tmpout;
    guidata(hObject,handles);
    uievscreen_plot(hObject,handles);
end
return;

% --- Executes on button press in ShowSpecBtn.
function ShowSpecBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ShowSpecBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowSpecBtn

fs=handles.fs;

[sm,sp,t,f]=evsmooth(handles.dat,fs,str2num(get(handles.SPTHBox,'String')));
axes(handles.SpecGramAxes);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);

handles.t=t;handles.f=f;
handles.sp=sp;
guidata(hObject,handles);
return;

% --- Executes on button press in QuitBtn.
function QuitBtn_Callback(hObject, eventdata, handles)
% hObject    handle to QuitBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.uievscreen);
return;

function SPTHBox_Callback(hObject, eventdata, handles)
% hObject    handle to SPTHBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SPTHBox as text
%        str2double(get(hObject,'String')) returns contents of SPTHBox as a double

[sm,sp,t,f]=evsmooth(dat,fs,str2num(get(handles.SPTHBox,'String')));
axes(handles.SpecGramAxes);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
return;

% --- Executes during object creation, after setting all properties.
function SPTHBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SPTHBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

% --- Executes on button press in ZoomBtn.
function ZoomBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom on;
return;

% --- Executes on button press in ZoomXBtn.
function ZoomXBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomXBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ZoomXBtn
zoom off;zoom xon;
return;

% --- Executes on button press in PanXBtn.
function PanXBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PanXBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pan xon;
return;

function FileNumBox_Callback(hObject, eventdata, handles)
% hObject    handle to FileNumBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileNumBox as text
%        str2double(get(hObject,'String')) returns contents of FileNumBox as a double
return;

% --- Executes during object creation, after setting all properties.
function FileNumBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileNumBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;
function varargout = evsonganaly(varargin)
% EVSONGANALY M-file for evsonganaly.fig
%      EVSONGANALY, by itself, creates a new EVSONGANALY or raises the existing
%      singleton*.
%
%      H = EVSONGANALY returns the handle to a new EVSONGANALY or the handle to
%      the existing singleton*.
%
%      EVSONGANALY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVSONGANALY.M with the given input arguments.
%
%      EVSONGANALY('Property','Value',...) creates a new EVSONGANALY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evsonganaly_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evsonganaly_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help evsonganaly

% Last Modified by GUIDE v2.5 19-Oct-2004 18:57:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @evsonganaly_OpeningFcn, ...
                   'gui_OutputFcn',  @evsonganaly_OutputFcn, ...
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
return;

% --- Executes just before evsonganaly is made visible.
function evsonganaly_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evsonganaly (see VARARGIN)

% Choose default command line output for evsonganaly
handles.output = hObject;

% get the structure array of input files
INPUTFILES=evloadfile;
if (length(INPUTFILES)==0)
    errordlg('No Input Files Found!');
    delete(handles.EVSONGANAL);
    return;
end
handles.INPUTFILES=INPUTFILES;
handles.NFILE=1;
handles.SPECTH=0.01;
handles.SEGTH=1.0e-4;
handles.MININT=5.0e-3;%in sec
handles.MINDUR=40.0e-3;%in sec
guidata(hObject, handles);handles=guidata(hObject);

%link the x axes
linkaxes([handles.SpecGramAxes,handles.LabelAxes,handles.SmoothAxes],'x');

% do first file
PlotDataFile(hObject,handles);
handles=guidata(hObject);

% setup some defaults for the GUI
set(handles.PrevFileBtn,'Value',get(handles.PrevFileBtn,'Min'));
set(handles.NextFileBtn,'Value',get(handles.NextFileBtn,'Min'));
set(handles.ResegBtn,'Value',get(handles.ResegBtn,'Min'));

%get initial zoom setting right
zoom xon;
set(handles.XZoomBtn,'Value',get(handles.XZoomBtn,'Max'));
set(handles.YZoomBtn,'Value',get(handles.YZoomBtn,'Min'));

uiwait(handles.EVSONGANAL);
return;

% --- Outputs from this function are returned to the command line.
function varargout = evsonganaly_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.EVSONGANAL);
return;

% --- Executes on button press in NextFileBtn.
function NextFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NextFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NextFileBtn

set(handles.NextFileBtn,'Value',get(handles.NextFileBtn,'Min'));
handles.NFILE=handles.NFILE+1;
if (handles.NFILE>length(handles.INPUTFILES))
    handles.NFILE=length(handles.INPUTFILES);
    errordlg('That was the last file there''s no going forward!');
end
guidata(hObject,handles);
PlotDataFile(hObject,handles);
handles=guidata(hObject);
return;

% --- Executes on button press in PrevFileBtn.
function PrevFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.PrevFileBtn,'Value',get(handles.PrevFileBtn,'Min'));
handles.NFILE=handles.NFILE-1;
if (handles.NFILE<1)
    handles.NFILE=1;
    errordlg('This is the first file there''s no going back!');
end
guidata(hObject,handles);
PlotDataFile(hObject,handles);
handles=guidata(hObject);
return;


% --- Executes on button press in ResegBtn.
function ResegBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ResegBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.ResegBtn,'Value',get(handles.ResegBtn,'Min'));
return;


% --- Executes on button press in QuitBtn.
function QuitBtn_Callback(hObject, eventdata, handles)
% hObject    handle to QuitBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.EVSONGANAL);
return;


% --- Executes on button press in XZoomBtn.
function XZoomBtn_Callback(hObject, eventdata, handles)
% hObject    handle to XZoomBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of XZoomBtn
val = get(hObject,'Value');
isYZoomOn = 0;
if (get(handles.YZoomBtn,'Value')==get(handles.YZoomBtn,'Max'))
        isYZoomOn = 1;
end
if (val==get(hObject,'Max'))
    if (isYZoomOn==1)
        zoom on;
    else
        zoom xon;
    end
else
    zoom off;
    if (isYZoomOn==1)
        zoom yon;
    end
end
return;


% --- Executes on button press in YZoomBtn.
function YZoomBtn_Callback(hObject, eventdata, handles)
% hObject    handle to YZoomBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of YZoomBtn

val = get(hObject,'Value');
isXZoomOn = 0;
if (get(handles.XZoomBtn,'Value')==get(handles.XZoomBtn,'Max'))
        isXZoomOn = 1;
end
if (val==get(hObject,'Max'))
    if (isXZoomOn==1)
        zoom on;
    else
        zoom yon;
    end
else
    zoom off;
    if (isXZoomOn==1)
        zoom yon;
    end
end
return;


% --- Executes on button press in LabelBtn.
function LabelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LabelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LabelBtn

onsets  = handles.ONSETS;
offsets = handles.OFFSETS;

% turn zoom off
zoom off;
set(handles.XZoomBtn,'Value',get(handles.XZoomBtn,'Min'));
set(handles.YZoomBtn,'Value',get(handles.YZoomBtn,'Min'));

CurXAxis = get(handles.SpecGramAxes,'XLim');
TxtLbl = handles.LABELTAGS;

pp = find(onsets>=CurXAxis(1));
if (length(pp)<1)
    CurLabelInd = 1;
else
    CurLabelInd = pp(1);
end
set(TxtLbl(CurLabelInd),'Color',[1,0,0]);

while (1)
    [XVal,YVal,KeyVal]=ginput(1);
    % CNTRL-Q or ESC to stop labeling
    if ((KeyVal==17)|(KeyVal==27))
        set(TxtLbl(CurLabelInd),'Color',[0,0,0]);
        break;
    end
    
    % left button press
    if ((KeyVal>=1)&(KeyVal<=3))
        axes(handles.SpecGramAxes);
        vv=axis;
        if (XVal>=vv(2))
            axis([vv(1:2)+diff(vv(1:2)),vv(3:4)]);
            vv=axis;
            XVal=vv(1);
        elseif (XVal<=vv(1))
            axis([vv(1:2)-diff(vv(1:2)),vv(3:4)]);
            vv=axis;
            XVal=vv(1);
        else
            if (KeyVal>1)
                dd = diff(vv(1:2));
                axis([XVal-dd,XVal+dd,vv(3:4)]);
            end
        end

        set(TxtLbl(CurLabelInd),'Color',[0,0,0]);
        pp = find(onsets>=XVal);
        if (length(pp)<1)
            CurLabelInd = 1;
        else
            CurLabelInd = pp(1);
        end
        set(TxtLbl(CurLabelInd),'Color',[1,0,0]);
    else
        %key press
        set(TxtLbl(CurLabelInd),'String',char(KeyVal));
        set(TxtLbl(CurLabelInd),'Color',[0,0,0]);
        if (CurLabelInd<length(TxtLbl))
            CurLabelInd = CurLabelInd + 1;
        else
            break;
        end
        set(TxtLbl(CurLabelInd),'Color',[1,0,0]);
    end
end

set(handles.LabelBtn,'Value',get(handles.LabelBtn,'Min'));
zoom xon;
set(handles.XZoomBtn,'Value',get(handles.XZoomBtn,'Max'));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DoNothing(src,evnt);
% does nothing;
return;
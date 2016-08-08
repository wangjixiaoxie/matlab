function varargout = lt_neural_sortAudioFiles(varargin)
%% ===
% varargin{1}='batch.keep';

% TO DO
% 1) deletion
% 2) stop when fid is -1


%% ==
% LT_NEURAL_SORTAUDIOFILES MATLAB code for lt_neural_sortAudioFiles.fig
%      LT_NEURAL_SORTAUDIOFILES, by itself, creates a new LT_NEURAL_SORTAUDIOFILES or raises the existing
%      singleton*.
%
%      H = LT_NEURAL_SORTAUDIOFILES returns the handle to a new LT_NEURAL_SORTAUDIOFILES or the handle to
%      the existing singleton*.
%
%      LT_NEURAL_SORTAUDIOFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LT_NEURAL_SORTAUDIOFILES.M with the given input arguments.
%
%      LT_NEURAL_SORTAUDIOFILES('Property','Value',...) creates a new LT_NEURAL_SORTAUDIOFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lt_neural_sortAudioFiles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lt_neural_sortAudioFiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lt_neural_sortAudioFiles

% Last Modified by GUIDE v2.5 04-Aug-2016 19:24:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lt_neural_sortAudioFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @lt_neural_sortAudioFiles_OutputFcn, ...
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


% --- Executes just before lt_neural_sortAudioFiles is made visible.
function lt_neural_sortAudioFiles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lt_neural_sortAudioFiles (see VARARGIN)

% Choose default command line output for lt_neural_sortAudioFiles
handles.output = hObject;

set(gcf,'toolbar','figure');

% ===== mksubfolder to delete songs
handles.subfolder='DELETED_SONGS';
if ~exist(handles.subfolder, 'dir')
    mkdir(handles.subfolder)
end

% ======= start reading batch
batchf=varargin{1};
fid=fopen(batchf);
handles.fid=fid;

% ====== start a new batch file for songs you want to keep
batchfnew=[batchf '.handKeep'];
fid2=fopen(batchfnew, 'w');
handles.fid2=fid2;
handles.batchfnew=batchfnew;

% ======== PLOT SPEC FOR SONG
filename=fgetl(fid);
axes(handles.axes1);
handles.plot1_filename=filename;
fn_PlotSpec(handles.plot1_filename);

filename=fgetl(fid);
axes(handles.axes8);
handles.plot8_filename=filename;
fn_PlotSpec(handles.plot8_filename);

filename=fgetl(fid);
axes(handles.axes9);
handles.plot9_filename=filename;
fn_PlotSpec(handles.plot9_filename);

filename=fgetl(fid);
axes(handles.axes10);
handles.plot10_filename=filename;
fn_PlotSpec(handles.plot10_filename);

filename=fgetl(fid);
axes(handles.axes11);
handles.plot11_filename=filename;
fn_PlotSpec(handles.plot11_filename);

filename=fgetl(fid);
axes(handles.axes12);
handles.plot12_filename=filename;
fn_PlotSpec(handles.plot12_filename);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lt_neural_sortAudioFiles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lt_neural_sortAudioFiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fid=handles.fid;
fid2=handles.fid2;
batchfnew=handles.batchfnew;

% ====== LT - RESET ALL FIGURES
% 1) pick out the ones you want to make new batch with
if get(handles.checkbox2, 'Value')==1
% --- write to new batch file    
    	fprintf(fid2,'%s\n',handles.plot1_filename);
        disp(['added ' handles.plot1_filename ' to batch ' batchfnew])
end

if get(handles.checkbox3, 'Value')==1
% --- write to new batch file    
    	fprintf(fid2,'%s\n',handles.plot8_filename);
        disp(['added ' handles.plot8_filename ' to batch ' batchfnew])
end

if get(handles.checkbox4, 'Value')==1
% --- write to new batch file    
    	fprintf(fid2,'%s\n',handles.plot9_filename);
        disp(['added ' handles.plot9_filename ' to batch ' batchfnew])
end

if get(handles.checkbox5, 'Value')==1
% --- write to new batch file    
    	fprintf(fid2,'%s\n',handles.plot10_filename);
        disp(['added ' handles.plot10_filename ' to batch ' batchfnew])
end

if get(handles.checkbox6, 'Value')==1
% --- write to new batch file    
    	fprintf(fid2,'%s\n',handles.plot11_filename);
        disp(['added ' handles.plot11_filename ' to batch ' batchfnew])
end

if get(handles.checkbox7, 'Value')==1
% --- write to new batch file    
    	fprintf(fid2,'%s\n',handles.plot12_filename);
        disp(['added ' handles.plot12_filename ' to batch ' batchfnew])
end

% === reset check boxes
        set(handles.checkbox2, 'Value', 0);
        set(handles.checkbox3, 'Value', 0);
        set(handles.checkbox4, 'Value', 0);
        set(handles.checkbox5, 'Value', 0);
        set(handles.checkbox6, 'Value', 0);
        set(handles.checkbox7, 'Value', 0);

        

% ========= DELETE SONGS
% (if radio button active, moves song to a subfolder)

if get(handles.radiobutton3, 'Value')==1
    % then move file to subfolder
    movefile(handles.plot1_filename, handles.subfolder);
end

if get(handles.radiobutton4, 'Value')==1
    % then move file to subfolder
    movefile(handles.plot8_filename, handles.subfolder);
end

if get(handles.radiobutton5, 'Value')==1
    % then move file to subfolder
    movefile(handles.plot9_filename, handles.subfolder);
end

if get(handles.radiobutton6, 'Value')==1
    % then move file to subfolder
    movefile(handles.plot10_filename, handles.subfolder);
end

if get(handles.radiobutton7, 'Value')==1
    % then move file to subfolder
    movefile(handles.plot11_filename, handles.subfolder);
end


if get(handles.radiobutton8, 'Value')==1
    % then move file to subfolder
    movefile(handles.plot12_filename, handles.subfolder);
end

% --------- RESET RADIO BUTTONS
set(handles.radiobutton3, 'Value', 0);
set(handles.radiobutton4, 'Value', 0);
set(handles.radiobutton5, 'Value', 0);
set(handles.radiobutton6, 'Value', 0);
set(handles.radiobutton7, 'Value', 0);
set(handles.radiobutton8, 'Value', 0);


 
% 2) replace all figs
% ======== PLOT SPEC FOR SONG
filename=fgetl(fid);
axes(handles.axes1);
handles.plot1_filename=filename;
fn_PlotSpec(handles.plot1_filename);

filename=fgetl(fid);
axes(handles.axes8);
handles.plot8_filename=filename;
fn_PlotSpec(handles.plot8_filename);

filename=fgetl(fid);
axes(handles.axes9);
handles.plot9_filename=filename;
fn_PlotSpec(handles.plot9_filename);

filename=fgetl(fid);
axes(handles.axes10);
handles.plot10_filename=filename;
fn_PlotSpec(handles.plot10_filename);

filename=fgetl(fid);
axes(handles.axes11);
handles.plot11_filename=filename;
fn_PlotSpec(handles.plot11_filename);

filename=fgetl(fid);
axes(handles.axes12);
handles.plot12_filename=filename;
fn_PlotSpec(handles.plot12_filename);



% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ===== ACTIOVATES ALL RADIO BUTTONS - ALL SONGS DELETED
set(handles.radiobutton3, 'Value', 1);
set(handles.radiobutton4, 'Value', 1);
set(handles.radiobutton5, 'Value', 1);
set(handles.radiobutton6, 'Value', 1);
set(handles.radiobutton7, 'Value', 1);
set(handles.radiobutton8, 'Value', 1);

% ==== then activates pushbutton 1 to reset all songs. [also ends this
% callback]
pushbutton1_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes8


% --- Executes during object creation, after setting all properties.
function axes9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes9


% --- Executes during object creation, after setting all properties.
function axes10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes10


% --- Executes during object creation, after setting all properties.
function axes11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes11


% --- Executes during object creation, after setting all properties.
function axes12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes12



% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7







% ----- Plots spec for this song
function fn_PlotSpec(fline)
% ===== if no more songs, then tell user, and quit program
if ~ischar(fline)
    plot([1 100], [1 10], 'r');
    msgbox('No More Songs! - can keep doing "Load Songs"');
else
    
    
    % =======================================================
    % LT - plot (Intan)
    [~, ~, frequency_parameters, board_adc_data]=pj_readIntanNoGui(fline);
    
    % --- NOTE ON THE OUTPUT VARIABLES
    
    
    % ---
    
    %
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title(fline);
    %     plot(board_adc_data, 'k');
    %     axis('tight');
    
    % === get params
    fs=frequency_parameters.board_adc_sample_rate;
    window = 0.016*fs; % make it 16ms, to match evtaf stuff
    nfft=2^nextpow2(window); % go to next pow 2 size
    olap=0.8;
    noverlap=olap*window;
    
    
    % ==== filter song
    songdat=board_adc_data;
    filter_type='hanningfir';
    F_low  = 500;
    F_high = 8000;
    songdat=bandpass(songdat,fs,F_low,F_high,filter_type);
    
    
    % === colect spectrogram
    [sp, f, t] = spectrogram(songdat, window, noverlap, nfft, fs);
    sp=abs(sp);
    
    % ========= cut off lower higher frequencies
    maxfreq=7500;
    inds=f<maxfreq;
    f=f(inds);
    sp=sp(inds, :);
    
    
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title(fline);
    
    
    % TAKE LOG of sp
    % first, convert any sp values of 0 to non-zero(to the lowest value present);
    % solves problem of taking log of 0
    pp=find(sp>0);
    mntmp = min(min(sp(pp)));
    pp=find(sp==0);
    sp(pp) = mntmp;
    
    % second, take log
    sptemp=log(sp);
    sptemp = sptemp - min(min(sptemp));
    sptemp = ((2^8) - 1)*(sptemp./max(max(sptemp)));
    
    % === black and red
    if (0)
        sptemp_vector=reshape(sptemp, numel(sptemp), 1);
        Xcenters=linspace(min(sptemp_vector), max(sptemp_vector), length(sptemp_vector)/200);
        [Ybinned, ~, ~]=lt_plot_histogram(sptemp_vector, Xcenters, 0, '','',1,'k');
        plot(Xcenters, Ybinned);
        
        [~, maxind]=max(Ybinned); % first peak in distrubtion of magnitudes
        Cmin=Xcenters(maxind);
        Cmax=max(sptemp_vector);
        
        % Plot
        colormap('hot')
        imagesc(t, f, sptemp, [Cmin+Cmin/10 Cmax]);
    end
    
    % Plot
    imagesc(t, f, sptemp);
    set(gca,'YDir','normal');
    title(fline);
end
%     axis([t(1) t(end) f(end) f(1)]);
%     axis([t(1) t(end) f(1) f(end)]);







% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called







% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8



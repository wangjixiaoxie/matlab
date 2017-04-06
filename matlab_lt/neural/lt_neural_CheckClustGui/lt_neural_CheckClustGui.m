function varargout = lt_neural_CheckClustGui(varargin)
%%
% RUN THIS IN NEURAL FOLDER.

% varargin{1} = 11; % channel num
% varargin{2} = 'Batch1438to1909'; % batch file
% varargin{3} = 14; % second chan. if empty [], then will not plot

%%
% LT_NEURAL_CHECKCLUSTGUI MATLAB code for lt_neural_CheckClustGui.fig
%      LT_NEURAL_CHECKCLUSTGUI, by itself, creates a new LT_NEURAL_CHECKCLUSTGUI or raises the existing
%      singleton*.
%
%      H = LT_NEURAL_CHECKCLUSTGUI returns the handle to a new LT_NEURAL_CHECKCLUSTGUI or the handle to
%      the existing singleton*.
%
%      LT_NEURAL_CHECKCLUSTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LT_NEURAL_CHECKCLUSTGUI.M with the given input arguments.
%
%      LT_NEURAL_CHECKCLUSTGUI('Property','Value',...) creates a new LT_NEURAL_CHECKCLUSTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lt_neural_CheckClustGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lt_neural_CheckClustGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lt_neural_CheckClustGui

% Last Modified by GUIDE v2.5 13-Mar-2017 19:38:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lt_neural_CheckClustGui_OpeningFcn, ...
                   'gui_OutputFcn',  @lt_neural_CheckClustGui_OutputFcn, ...
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


% --- Executes just before lt_neural_CheckClustGui is made visible.
function lt_neural_CheckClustGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lt_neural_CheckClustGui (see VARARGIN)

% =============== PLOT
% hold on;
% x = rand(20);
% y = rand(20);
% hplot = plot(x, y, 'ok');

%% --------------------------------------
if (1)
channel_board = varargin{1};
batchf = varargin{2};
plotRawSound = 1;
plotcols={'m', 'r','c', 'b', 'g'};

    PlotSecondChan=0;
if ~isempty(varargin{3})
    PlotSecondChan=1;
    SecondChan = varargin{3};
end

% --- go to data folder
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)

% ---- load, go back up, and plot.
% 1) EXTRACT AND CONCAT ALL SONG FILES and neural files
load('MetaDat.mat'); % contains file names
% --- load concat neural and spikes
neural_cat=load('data.mat');
spikes_cat=load('times_data.mat');


AllSongs_old=[];
AllNeural_old=[];
AllLabels=[];
AllOnsets=[];
AllOffsets=[];
cumulative_filedur=0; % keeps track as concatenating

AllNeural_SecondChan = [];

cd ..

cumulativetime_ms = 0;
% --- load audio and old neural (from actual files)
mm=1;
    hsplots = [];
    % -- load original sound and neural
    [amplifier_data,~,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(mm).filename);
    
    indsamp=find([amplifier_channels.chip_channel]==channel_board);
    neurdat = amplifier_data(indsamp, :);
    
    % -- load labels, onsets, offsets
    tmp=load([metaDat(mm).filename '.not.mat']);
    AllLabels=[AllLabels tmp.labels];
    
    tt = [1:length(board_adc_data(1,:))]/frequency_parameters.amplifier_sample_rate;
    
    datdur_ms = 1000*length(board_adc_data(1,:))/frequency_parameters.amplifier_sample_rate;
    
    axes(handles.axes1)
    hold on;
    
    % - a. song raw
hsplot = lt_subplot(2,1,1); hold on;
plot(tt, board_adc_data(1,:));
        hsplots=[hsplots hsplot];
    
    % onsets and labels
    for i=1:length(tmp.onsets)
        line([tmp.onsets(i) tmp.offsets(i)]/1000, [1.65 1.65], 'LineWidth', 2, 'Color', 'r');
        lt_plot_text(tmp.onsets(i)/1000, 1.5, tmp.labels(i), 'r')
    end
    
    
    if PlotSecondChan==1
        ind =find([amplifier_channels.chip_channel]==SecondChan);
hsplot = lt_subplot(2,1,2); hold on;
[datfilt] =lt_neural_filter(amplifier_data(ind, :), frequency_parameters.amplifier_sample_rate);
        plot(tt, datfilt, 'k');
        hsplots=[hsplots hsplot];
        title(['second chan (' num2str(SecondChan) ')']);
        
    end
    
    
    % - d. neural (cat)
    axes(handles.axes2)
    hold on;
% hsplot = lt_subplot(2,1,2); hold on;
[datfilt] =lt_neural_filter(neurdat, frequency_parameters.amplifier_sample_rate);
    plot(tt, datfilt, 'k');
    hsplots=[hsplots hsplot];

    
    % overlay spike times (cat)
    numclust=max(spikes_cat.cluster_class(:,1));
    for i=1:numclust
        inds=find([spikes_cat.cluster_class(:,1)]==i & ...
            spikes_cat.cluster_class(:,2)>cumulativetime_ms ...
            & spikes_cat.cluster_class(:,2)<=cumulativetime_ms+datdur_ms);
        
        for iii=1:length(inds)
            ii=inds(iii);
            spktime=spikes_cat.cluster_class(ii, 2); % in ms
            spktime=spktime-cumulativetime_ms;
            line([spktime spktime]/1000, [0 40]-i*20, 'Color', plotcols{i}, 'LineWidth',2);
        end
    end
    
    
    % --- link all
    linkaxes(hsplots, 'x');
    
    %         if PlotSecondChan==1
    %             ind=find([amplifier_channels.chip_channel]==SecondChan);
    %             AllNeural_SecondChan=[AllNeural_SecondChan amplifier_data(ind, :)];
    %         end
    cumulativetime_ms = cumulativetime_ms + datdur_ms;
    
%     lt_subtitle(['song #' num2str(mm) '/' num2str(length(metaDat)) '(' metaDat(mm).filename ')'])
% --------------------------
end

%% --------------------------------------
if (0)
    channel_board = varargin{1};
batchf = varargin{2};
plotRawSound = 1;
plotcols={'m', 'r','c', 'b', 'g'};

    PlotSecondChan=0;
if ~isempty(varargin{3})
    PlotSecondChan=1;
    SecondChan = varargin{3};
end

% --- go to data folder
datdir=['Chan' num2str(channel_board) 'amp-' batchf];
cd(datdir)

% ---- load, go back up, and plot.
% 1) EXTRACT AND CONCAT ALL SONG FILES and neural files
load('MetaDat.mat'); % contains file names
% --- load concat neural and spikes
neural_cat=load('data.mat');
spikes_cat=load('times_data.mat');


AllSongs_old=[];
AllNeural_old=[];
AllLabels=[];
AllOnsets=[];
AllOffsets=[];
cumulative_filedur=0; % keeps track as concatenating

AllNeural_SecondChan = [];

cd ..

cumulativetime_ms = 0;
% --- load audio and old neural (from actual files)
mm=1;
    hsplots = [];
    % -- load original sound and neural
    [amplifier_data,~,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(metaDat(mm).filename);
    
    indsamp=find([amplifier_channels.chip_channel]==channel_board);
    neurdat = amplifier_data(indsamp, :);
    
    % -- load labels, onsets, offsets
    tmp=load([metaDat(mm).filename '.not.mat']);
    AllLabels=[AllLabels tmp.labels];
    
    tt = [1:length(board_adc_data(1,:))]/frequency_parameters.amplifier_sample_rate;
    
    datdur_ms = 1000*length(board_adc_data(1,:))/frequency_parameters.amplifier_sample_rate;
    
    
    
    % - d. neural (cat)
    hold on;
    subplot(3,1,3); hold on;
[datfilt] =lt_neural_filter(neurdat, frequency_parameters.amplifier_sample_rate);
    plot(tt, datfilt, 'k');

%     title('neural, concated file (extracted spikes)');
    
    % overlay spike times (cat)
    numclust=max(spikes_cat.cluster_class(:,1));
    for i=1:numclust
        inds=find([spikes_cat.cluster_class(:,1)]==i & ...
            spikes_cat.cluster_class(:,2)>cumulativetime_ms ...
            & spikes_cat.cluster_class(:,2)<=cumulativetime_ms+datdur_ms);
        
        for iii=1:length(inds)
            ii=inds(iii);
            spktime=spikes_cat.cluster_class(ii, 2); % in ms
            spktime=spktime-cumulativetime_ms;
            line([spktime spktime]/1000, [0 40]-i*20, 'Color', plotcols{i}, 'LineWidth',2);
        end
    end
    
    if PlotSecondChan==1
        ind =find([amplifier_channels.chip_channel]==SecondChan);
hsplot = lt_subplot(3,1,3); hold on;
[datfilt] =lt_neural_filter(amplifier_data(ind, :), frequency_parameters.amplifier_sample_rate);
        plot(tt, datfilt, 'k');
        hsplots=[hsplots hsplot];
        title(['second chan (' num2str(SecondChan) ')']);
        
    end
    
    % --- link all
%     linkaxes(hsplots, 'x');
%     
%     %         if PlotSecondChan==1
%     %             ind=find([amplifier_channels.chip_channel]==SecondChan);
%     %             AllNeural_SecondChan=[AllNeural_SecondChan amplifier_data(ind, :)];
%     %         end
%     cumulativetime_ms = cumulativetime_ms + datdur_ms;
%     
%     lt_subtitle(['song #' num2str(mm) '/' num2str(length(metaDat)) '(' metaDat(mm).filename ')'])
% --------------------------
end

% Choose default command line output for lt_neural_CheckClustGui
% hObject = gca;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lt_neural_CheckClustGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lt_neural_CheckClustGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% 
% % --- Executes on mouse press over axes background.
% function axes1_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to axes1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% set(gcf, 'Units', 'Normalized')
% disp('WAITING!!!')
% waitforbuttonpress
% finalrect = rbbox;
% 
% annotation('rectangle', finalrect);
% 
% %===== convert from figure units to units in the figure
% set(gca, 'Units', 'normalized')
% axes_pos = get(gca, 'Position'); % position of axes rel to figure (norm)
% 
% Xlim = get(gca, 'XLim'); % axes limits (actual units) 
% Ylim = get(gca, 'Ylim');
% xrange = Xlim(2)-Xlim(1);
% 
% % -- actual unit for left limit
% tmp = (finalrect(1) - axes_pos(1))/axes_pos(3);
% X1 = Xlim(1)+tmp*xrange
% 
% plot(X1, 0.3, 'og');
% 
% % -- actual unit for right limit
% tmp = (finalrect(1)+finalrect(3) - axes_pos(1))/axes_pos(3);
% X2 = Xlim(1) + tmp*xrange
% 
% plot(X2, 0.3, 'or');
% 
% 
% 
% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf, 'Units', 'Normalized')
disp('WAITING!!!')
waitforbuttonpress
finalrect = rbbox;

annotation('rectangle', finalrect);

%===== convert from figure units to units in the figure
set(gca, 'Units', 'normalized')
axes_pos = get(gca, 'Position'); % position of axes rel to figure (norm)

Xlim = get(gca, 'XLim'); % axes limits (actual units) 
Ylim = get(gca, 'Ylim');
xrange = Xlim(2)-Xlim(1);

% -- actual unit for left limit
tmp = (finalrect(1) - axes_pos(1))/axes_pos(3);
X1 = Xlim(1)+tmp*xrange

plot(X1, 0.3, 'og');

% -- actual unit for right limit
tmp = (finalrect(1)+finalrect(3) - axes_pos(1))/axes_pos(3);
X2 = Xlim(1) + tmp*xrange

plot(X2, 0.3, 'or');

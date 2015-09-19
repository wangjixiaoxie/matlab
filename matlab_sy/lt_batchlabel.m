%% LT 10/28 - changed line 103 to SegmentNotesJC, since that is used in evsonganaly (instead of evsegment)
% also line 117 from feature_vect to lt_feature_vect (conv, not smooth,
% used)
% also line 108, multiplied on and offsets by 1000, since subsequent analysis, and saving of notmat fiel, requires to be in ms format.
% also, changed the way songs segmented, so that it exactly matches
% evsonganaly: 1) ReadDataFile line 63; 2) SmoothData line 115; 3)
% SegmentNotesJC line 120. (previous version gave times that were 8ms ahead
% of evsonganalys.) - Note, cause of that lag is lack of filtfilt.


function lt_batchlabel(songfile, templates, template_cov, template_mean, index)
% Labels a song and saves the labels in a .not.mat file. 
% changed clint's version to not save all the files since
% disk space is becoming scarce, also do not need to run
% bactch uisonganal in order to use this version

global labelIndex
global filetype
global label_path
global id_threshold

global GL_THRESHOLD
global GL_MIN_INT
global GL_MIN_DUR
global GL_SM_WIN
global REUSENOTMAT


disp(['Labeling ' songfile]);

baresongfile=songfile;
if (~(songfile(1)=='/'))&(~strcmp(label_path,pwd))
     songfile=[label_path,songfile];
end

if (songfile(end-4:end)=='.filt')
     filtfile = songfile;
     songfile = songfile(1:end-5);
     if ~strcmp(filetype,'filt')
         disp('Filetype .filt found.')
     end
     filetype='filt';
else
     filtfile = [songfile, '.filt'];         
     if (songfile(end-4:end)=='.cbin')
           if ~(strcmp(filetype,'obs0r')|strcmp(filetype,'obs1r'))
               disp('Filetype .cbin found.')
           end
           filetype='obs0r';
     else
           if ~(strcmp(filetype,'w'))
               disp('Assuming wave file format (file does not end in .cbin or .filt)')
           end
           filetype='ebin0r';
     end
end
notefile = [songfile, '.not.mat'];


% If filtfile exists, read it, otherwise filter rawsong and save
% if exist(filtfile)
%     [filtsong, fs]=read_filt(filtfile);
% else
    if exist(songfile)
%         [song,fs]=evsoundin(label_path, baresongfile, filetype); 
    [song,fs,DOFILT,ext]=ReadDataFile([label_path baresongfile],'0'); % cbin, 0
        
% 	disp(['fs=',num2str(fs)]);
% 	disp(filetype);
%         filtsong=bandpass(song,fs,300,8000);
        %write_filt(filtfile, filtsong, fs);
        %disp(['    Filtering ',songfile])
    else
        disp(['ERROR:  Neither ',songfile,' nor ',filtfile,' can be found.'])
        return;
    end
% end

% if notefile exists use its segmentation boundries
if (exist(notefile)&(REUSENOTMAT))
    load(notefile);
    if length(onsets)==0
        disp(['No syllables in file ',notefile])
        return
    end

    % DON"T save old notefile 
    %eval(['save ', notefile, '.old', ...
    %   ' Fs',' onsets',' offsets',' labels',...
    %   ' threshold',' min_int',' min_dur',' sm_win']);

    on  = fix(onsets.*fs/1000);
    off = fix(offsets.*fs/1000);
%     if off(end) > length(filtsong)
%         off(end) = length(filtsong);
%     end
    if off(end) > length(song)
        off(end) = length(song);
    end
  
else
    %disp(['The notefile -- ',notefile,' -- does not exist!']);
    %return;
    
    Fs = fs;
    threshold = GL_THRESHOLD;
    min_int = GL_MIN_INT;
    min_dur = GL_MIN_DUR;
    sm_win  = GL_SM_WIN;
    
    %squared_song = filtsong.^2;
    %len=round(Fs*sm_win/1000.0);hhh=ones(1,len)/len;
    %smooth=conv(hhh,squared_song);
    %doffset=round((length(smooth)-length(filtsong))/2);
    %smooth=smooth(1+doffset:length(filtsong)+doffset);
    
%     [smooth]=evsmooth(filtsong,Fs,0.01);
    [smooth,sp,t,f]=SmoothData(song,Fs,1,'hanningfirff');
    
%     disp(['th=',num2str(threshold)]);
%     [onsets, offsets] = SegmentNotes(smooth,fs,min_int,min_dur,threshold);
    [onsets, offsets] = SegmentNotesJC(smooth,fs,min_int,min_dur,threshold);
    onsets=onsets*1000;
    offsets=offsets*1000;
    on  = floor(onsets.*fs/1000.0)+1;
    off = floor(offsets.*fs/1000.0)+1;

    %eval(['save ',notefile,' Fs onsets offsets labels threshold ',...
    %      'min_int min_dur sm_win']);
end

% calculate feature vectors
featvect = [];
% for i = 1:size(on,1)
%     start = on(i);
%     stop = off(i);
%     if ((start<length(filtsong))&(stop<=length(filtsong)))
% 	    featvect = [featvect;lt_feature_vect(filtsong(start:stop),32000)];
%     end
% end

for i = 1:size(on,1)
    start = on(i);
    stop = off(i);
    if ((start<length(song))&(stop<=length(song)))
	    featvect = [featvect;lt_feature_vect(song(start:stop),32000)];
    end
end

nlabels=zeros(size(featvect,1), 1); % labels as numerical index to templates array

% size of featvect, number of unknowns
sfv=size(featvect,1);
% number of elements in index, number of syllable labels 
maxi=max(index);

%% calculate mahalonobis distance to each group
%% this is the core code from the matlab function 'classify'. that function isn't used so that some information
%% can be saved for the confidence measure.
%d = zeros(sfv,maxi);
%for k = 1:maxi
%    templatek = templates(find(index == k),:);
%    d(:,k) = mahal(featvect,templatek);
%end

if (length(on)>0)
for k = 1:maxi
    kmean=template_mean(k,:);
    kdiff=(featvect-kmean(ones(size(featvect,1),1),:));
    kdiff(:,6)=kdiff(:,6)+100*(kdiff(:,6)<-50)-100*(kdiff(:,6)>50);
    kdiff(:,7)=kdiff(:,7)+100*(kdiff(:,7)<-50)-100*(kdiff(:,7)>50);
    kinvcov=inv(template_cov{k});
    d(:,k) = sqrt(sum((kdiff*kinvcov).*kdiff,2));
end

% class is the group to which the mahal distance is smallest
[mindist, class] = min(d');    % mindist == the smallest mahal distance -- use this to find syllables
class = class';                %          which are far from even the closest template

labels=[];

% Turn number into letter label
for i = 1:size(featvect,1)    
    labels=[labels; labelIndex(class(i,1))];
end

labels=labels';   % to match format of uisonganal

labels(find(mindist>id_threshold*size(featvect,2)))='-';  % Mark these as unknowns, otherwise mark everything with the best match.
else 
	labels = [];
end

% save results
eval(['save ', notefile,...
       ' Fs',' onsets',' offsets',' labels',...
       ' threshold',' min_int',' min_dur',' sm_win']);


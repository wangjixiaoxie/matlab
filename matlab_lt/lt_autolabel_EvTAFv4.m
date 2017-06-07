function [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat)
%% LT 5/18/15 - Autolabel song files by using simulation of evtafv4. given template files (in .evconfig format), assigns any detects as the desired syllable
% e.g. inputs:
% batch = 'batch.labeled.all';
% config= '/bluejay4/lucas/birds/rd23gr89/042815_CtxtDepPitch2_TargB2__day3/config.evconfig2';
% syl is structure telling how to fill in syls.
% syl.targ='b' label trigger
% syl.pre='cc' string preceding trigger
% syl.post='dd' string post trigger
% NoteNum= which evtafv4 note whose triggers to use? (0, 1, ...)
% Segmentation params:
% ampThresh, min_dur, min_int (as in evsonganaly)
% overwrite_notmat=1 (will overwrite old not mat. i.e. will save copy of
% notmat in backup folder (only for songs in this batch) and create new
% notmat.
% overwrite_notmat=0 (if notmat exists, will add onto it. i.e. keeps old
% labels, unless that label conflicts with new label). also saves a
% backup copy of old notmat

% OUTPUTS:
% [fnames, sylnum, vlsorfn, vlsorind] are used for 
% lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind) after checking wav file using evsonganaly - see lt_autolabel_MasterScript

% NOTE: .wav fiel it makes using only single syl, not context

% Run this in the folder containing all the songs.


%% Default params

CHANSPEC='obs0'; % for opening labview files
sm_win=2; % used for saving notmat

SylBoundsAdd=16; % in ms, how much time to remove from onset and add to offset in permissive window for labeling? I make it 0.016s (2 fft bins)
% this is useful for the time when a syl is fragmented (e.g. problem
% with thresholding), or trig is outside onset or offset by a little.

% How many syls in string used for labeling?
NumLabelsPre=length(syl.pre);
NumLabelsPost=length(syl.post);
LabelStr=[syl.pre syl.targ syl.post];


%% make folder to save backups of old notmats
tstamp=lt_get_timestamp(0);

OldNotMatdir=['OldNotMat_AL_moved_' tstamp];
mkdir(OldNotMatdir);


%% Run simulation of evtaf to get detects
disp('Running evtaf sim...')
AllData=EvTAFv4Sim_LT_temp(batch, config,'obs0');

disp(['Found and simulated evtaf for ' num2str(length(AllData)) ' songs today']);


%% Go through every song in batch and look at ttimes extracted above, and replace those segmented notes with syl letter, if it contains a trigger


fid=fopen(batch,'r');
count = 1;

while (1)
    fn=fgetl(fid);
    disp(fn)
    
    if (~ischar(fn));
        break;
    end;
    
    if (~exist(fn,'file'))
        continue;
    end
    
    
    % === Extract ttimes for this song
    % assume count num corresponds to index in AllData, but will also check to
    % make sure names of files match
    if strcmp(AllData(count).filename, fn)==0;
        disp('PROBLEM - filename of extracted data does not match batch file...')
        keyboard
    end
    
    % -- If this song has no triggers, still make notmat file, but just
    % won't label.
    if ~isempty(AllData(count).data_OfflineAllDetects);
        
        % -- get trigger times
        ttimes=AllData(count).data_OfflineAllDetects(:,1);
        
        % -- Only keep ttimes for desired note
        notenums=AllData(count).data_OfflineAllDetects(:,2);  
        ttimes=ttimes(notenums==NoteNum);
    else
        ttimes = [];
    end
    

    % === Look in song for syllable that corresponds to thiose ttimes
    
    % -- Make, or extract .not.mat info, depending on whetehr want to overwrite
    % old notmat (and whether it exists)
    % - does a notmat already exist?
    if (~exist([fn,'.not.mat'],'file')); % does not already exist
        
        % then create new notmat
        [datFull,fs]=evsoundin('',fn,CHANSPEC);
        sm=SmoothData(datFull,fs,1,'hanningfirff');
        sm(1)=0.0; sm(end)=0.0;
        [ons,offs]=SegmentNotesJC(sm,fs,min_int,min_dur,ampThresh);
        onsets=ons*1e3;offsets=offs*1e3;
        labels = char(ones([1,length(onsets)])*fix('-'));
    else
        % then notmat exists - first save backup of old notmat
        eval(['!cp ' fn '.not.mat ' OldNotMatdir]);
        
        % second, decide whether to replace, or simply append
        if overwrite_notmat==1;
            % delete old notmat
            eval(['!rm ' fn '.not.mat'])
            
            % then create new notmat
            [datFull,fs]=evsoundin('',fn,CHANSPEC);
            sm=SmoothData(datFull,fs,1,'hanningfirff');
            sm(1)=0.0;sm(end)=0.0;
            [ons,offs]=SegmentNotesJC(sm,fs,min_int,min_dur,ampThresh);
            onsets=ons*1e3;offsets=offs*1e3;
            labels = char(ones([1,length(onsets)])*fix('-'));
            
        else
            % load old notmat -
            load([fn,'.not.mat']);
            fs = Fs;
        end
    end
    
    
    % -- Now have notmat info.  Stick syllables in there
    for ii = 1:length(ttimes)
        pp=find(((onsets-SylBoundsAdd)<=ttimes(ii)) & ((offsets+SylBoundsAdd)>=ttimes(ii))); % notes coinciding with ttime
        
        if (length(pp)>0); % multiple notes...
            %                 labels(pp(1))=sylLabel; % greater than 1 is almost certainly multiple trigs, unless SylBoundsAdd was set larger than smallest gap dur, take 1st syl with trig.
            try % try to stick flanking sequence in. if fails, then fits as much as possible on either side.
                labels(pp(1)-NumLabelsPre:pp(1)+NumLabelsPost)=LabelStr; % greater than 1 is almost certainly multiple trigs, unless SylBoundsAdd was set larger than smallest gap dur, take 1st syl with trig.
            catch err
                disp(['error: song ' num2str(ii) ' could not fit entire label string, will fill with what I can']);
                
                % find out how much of label can stick in there.
                if pp(1)-1<NumLabelsPre;
                    Npre=pp(1)-1;
                else
                    Npre=NumLabelsPre;
                end
                if length(onsets)-pp(1)<NumLabelsPost;
                    Npost=length(onsets)-pp(1);
                else
                    Npost=NumLabelsPost;
                end
                
                % fill in labels as can
                labels(pp(1)-Npre:pp(1)+Npost)=LabelStr(NumLabelsPre-Npre+1:end-(NumLabelsPost-Npost));
            end
        end
    end
    
    
    % == Save notmat and update counter
    fname=fn;
    Fs = fs;
    threshold = ampThresh;
    cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
        'offsets onsets sm_win threshold'];
	eval(cmd);

    clear sm

    count = count + 1;
end


%%  == Save new version of batch file with "autolabeled" in name
batch_new=[batch '.AutoLabeled'];
eval(['!cp ' batch ' ' batch_new])


%% ======================= TO CHECK ACCURACY - isolates all syls labeled and saves to .wav - open and check by eye.
[fnames, sylnum]=lt_jc_chcklbl(batch, syl.targ, 0.025,0.025,'','','');
[vlsorfn vlsorind]=jc_vlsorfn(batch, syl.targ,'','');

% Troubleshooting, enter syl name you want to find song of.
% [fnames, sylnum]=lt_jc_chcklbl(batch,'x', 0.025,0.025,'','','');
% [vlsorfn vlsorind]=jc_vlsorfn(batch,'x','','');

disp('DONE!');




disp('Done autolabeling!')

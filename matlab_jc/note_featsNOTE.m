function [notefeats, mean_syllfeats, std_syllfeats, labelIndex]=note_featsNOTE(batch, filetype);

%%modified from autolabel_cara
%%reads in a batch file of .wavs
%and calculates features for all syllables
%saves matrix of notefeats, all features for all syllables
%%saves matrixes of mean/std of 8 features of each syllable in labelIndex
%%currently 8 features from feature_vect_test (see below)
     
         
sylls=[];
labelIndex=[];
notefeats=[];
index=[];
id=[];
templates = [];
all_labels=[];
counts = zeros(size(labelIndex));
mean_syllfeats_current = zeros([size(labelIndex,1),(size(notefeats,2))]);
mean_syllfeats = zeros([size(labelIndex,1),(size(notefeats,2))]);
std_syllfeats = zeros([size(labelIndex,1),(size(notefeats,2))]);

%open batch_file of .wavs
fid=fopen(batch,'r');
%get name of metafile containing notefile names
while (1)
    fn = fgetl(fid);
    if (~ischar(fn))
        break;
    end
notefile=[fn, '.not.mat'];
filtfile=[fn, '.filt'];
% Read filt file if it exists, otherwise read and filter raw song, saving a copy
    if exist(filtfile)
       [filtsong, fs]=read_filt(filtfile);
    else
        disp(['   Filtering ',fn])
        [song,fs]=soundin('', fn, filetype);
        filtsong=bandpass(song,fs,300,8000);
        write_filt(filtfile, filtsong, fs);
    end
%Read in .not.mat variables   
    load(notefile);
    %make sure labels is a column vector
    [x,y]=size(labels);
    if y > x
        labels=labels';
    end    
    disp(['Processing ' fn ' . . .']);     
    % convert onsets and offsets from ms to samples
    on = fix(onsets.*fs/1000);
    off = fix(offsets.*fs/1000);
  
    % read in labels and make index of labels found
    all_labels=[all_labels; labels];
    % for each syllable in the song file, get timing info and label info
    for ii = 1:size(labels,1)   
        if labels(ii)=='b'
        %add time of current syllable to 'sylls' vector
        sylls = [sylls; {filtsong(on(ii):off(ii))}];
        %add syllable label to labelindex if it isn't there already
        if isempty(labelIndex) | isempty(find(labelIndex == labels(ii)))
            labelIndex = [labelIndex; labels(ii)];
        end
        %add syllable label to index of all syllables (sylls), in numerical form
        index = [index; find(labelIndex == labels(ii))];
        id = [id; {fn},{onsets(ii)}];
        end
    end
end
    %Calculate features for every syllable present using feature_vect_test
    disp(['Calculating features for all  ', num2str(size(sylls,1)), '  syllables, go get some coffee']);
    for ii = 1:size(sylls, 1);  
        current_feature=feature_vect_test(sylls{ii},fs);    
        notefeats = [notefeats; current_feature];
    end
        % Count total number of examples of syllables from all files and
        % report results
    for ii = 1:size(labelIndex, 1);
        counts(ii) = size(find(index == ii),1);
        disp(['Found ' int2str(counts(ii)) ' examples of ' labelIndex(ii)]);
        syllable_index=(find(index == ii));
        %notefeats(a_index,5)
        %find indexes for each label, and calculate 8 features for each
        %syllable separately
        for xx=1:size(notefeats,2);
            current_mean=mean(notefeats(syllable_index,xx));
            current_std=std(notefeats(syllable_index,xx));
            disp(['F', num2str(xx)  'syllable ', labelIndex(ii), '  =  ', num2str(current_mean) ' +/-  ', num2str(current_std)]) ,
            mean_syllfeats(ii,xx)=current_mean;
            std_syllfeats(ii,xx)=current_std;
        end
    end

    
    %display means for all syllables, 8 features total
    disp(['Means for 8 features in   ', batch, ' are: ']);
    for mm=1:size(notefeats,2)
        disp(['all syllables   ',num2str(mm), '  =  ',num2str(mean(notefeats(:,mm))),'     +/-  ', num2str(std(notefeats(:,mm))),'     N =  ', num2str(length(notefeats(:,mm)))]);    

    end
end

       
    

% (see feature_vect_test)
   % 1: Entropy of the spectral density 
   % range: 0 to 100; 100 = white noise; 0 = pure note
   %featvect(1)=-100*sum(sd.*log2(sd))/log2(length(sd));
   % 2: Duration of Syllable (in ms);
   %featvect(2)=10*log2(Dur*1000);
   % 3: Entropy of the loudness vs. time 
   %(range: 0 to 100; 100 = flat amplitude; 0 = one very sharp peak)
   %featvect(3)=-100*sum(loud.*log2(loud))/log2(length(loud));
   % 4: Amplitude Slope;
   %featvect(4)=50+50*(sum1-sum2)/(sum1+sum2);  
   % 5: spectrotemp entropy; diff between sweep and stack
   %featvect(5)=100*spec_ent;
   % 6: Time to half peak loudness
   %featvect(6)=100*halftime;
   % 7: The Mean Frequency; dependent on loudness?
   %featvect(7)=sum(sd.*sdf)/100;
   % 8: The Frequency Slope
   %featvect(8)=50+50*(f1-f2)/(f1+f2);



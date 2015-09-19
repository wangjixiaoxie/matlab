function autolabel()

% AUTOLABEL()
%
% Automatically labels the syllables of a birdsong. This program prompts the user for a batch file with the names
% of hand labeled songs as well a batch file with the names of the songs to be labeled. The hand labeled songs should 
% have at least 30 examples of every syllable type in a .not.mat file of the same name. If a syllable type has less 
% than 20 examples it will be ignored and this program will not attempt to give any syllables that label. That means
% that all syllables of that type will be mislabeled or left blank. This program also ignores any syllables marked 
% '0' or '-'.  Marking ambiguous or unclear syllables as '0' in your labelling will help the program to find distinct
% clusters.
%
% Syllables that cannot be confidently labeled are left blank (given the label '-'). After labeling is done AUTOLABEL 
% will optionally run LABELBLANKS, an interface for quickly hand labeling syllables with the label '-'.
% LABELBLANKS can be run seperately and will prompt the user for a batch file of songs with '-' labeled syllables.
% That will be the same batch file of songs given to AUTOLABEL.
%
% This program WILL make mistakes in addition to leaving syllables blank, so don't trust it too much. 
%

global labelIndex %array of the letters of labels in use

global filetype
global id_threshold

%filetype = get_filetype;
filetype='obs0r';
filetype='w';
trainup=-1;
while (trainup<0)
    disp('Would you like to:');
    disp('     (1)   <L>oad a template file');
    disp('  or (2)   <T>rain on a set of labeled files?');
    tmp=input('                ','s');
    if(length(tmp)==0)
        tmp='qq';
        disp(' Please select 1 or 2. ');
    end
    if (tmp(1)=='L') | (tmp(1)=='l') | (tmp=='1')
        trainup=0;
    end
    if (tmp(1)=='T') | (tmp(1)=='t') | (tmp=='2')
        trainup=1;
    end
end

% Get batchfiles. training_fid is labeled batchfile, songs_fid is unlabeled batchfile
if trainup
    [training_fid] = get_training_batchfile;    
else
    disp('Select the template file to be used.')
    uiload;
    if (~exist('templates'))|(~exist('labelIndex'))|(~exist('index'))
        disp('ERROR: The template file does not contain the needed variables.')
        return
    end
    if (~exist('template_cov'))|(~exist('template_mean'))
        disp('   OLD TEMPLATE FILE -- calculating covariance matrices...')
        [template_cov,template_mean] = get_covariances(templates,index);
        disp('   saving new version of template file as "template.mat".')
        save 'template.mat' labelIndex index templates template_cov template_mean;
    end
end    
[songs_fid] = get_labelling_batchfile;
pause(.1); % to let popup menu clear

disp('What threshold would you like to use, to determine which syllables are unidentifiable?')
tmp=input('Higher thresholds give fewer unlabelled syllables (default=1.5) ','s');
if length(tmp)==0
    tmp='1.5';
end
id_threshold=str2num(tmp);
if ~exist('id_threshold')|(id_threshold==NaN)|(id_threshold<=0)
    id_threshold=1.5;
end
disp(['     using a threshold of ',num2str(id_threshold)])

% Read info from labeled files. labelIndex is an array of the labels that will be used. Their position in 
% the array is their numerical equivalent. 
if trainup
    [templateSylls, index, labelIndex, temp_id] = read_training(training_fid);
end

% Makes the template and template structure, then save
if trainup
    templates = make_templates(templateSylls);
    [template_cov,template_mean] = get_covariances(templates,index);
    strctname = input('What would you like to name this template structure? ' );  
    eval([char(strctname) '.labelIndex = labelIndex;']); 
    eval([char(strctname) '.index = index;']); 
    eval([char(strctname) '.templates = templates;']);
    eval([char(strctname) '.temp_id = temp_id;']); 
    eval([char(strctname) '.template_cov = template_cov;']); 
    eval([char(strctname) '.template_mean = template_mean;']);
    sv_name = [char(strctname) '.template_strct.mat'];
    save(char(sv_name), char(strctname));
end                     

tmp=[];
if ~trainup
    tmp=input('Would you like to evaluate the template? ','s');
    if length(tmp)==0
        tmp='No';
    end
end

lblornot  = input('Do you want to evaluate the templates now and continue with the autolabeling process? [0 or 1]');

if lblornot
        
	
	% Evaluate the templates -- how tight is each syllable definition, and how
	% distinct are the different syllables.
	if trainup | tmp(1)=='Y' | tmp(1)=='y'
        evaluate(templates,template_cov,template_mean,index);
	end
	
	% Label each song in songs batchfile
	songfile = fscanf(songs_fid,'%s',1);
	while ~isempty(songfile)
         batchlabel(songfile, templates,template_cov,template_mean, index); 
         songfile = fscanf(songs_fid,'%s',1);
	end
	disp('End of songfile.')
	
	% Display some stats about what the program has done 
	countblanks(songs_fid);
            
	% Optionally run labelblanks.
	%    disp('Would you like to label the syllables that were left blank now? (y/n) [y]');
	%    tmp = input(' ', 's');
	%    if strcmp(tmp, 'y') | strcmp(tmp, 'Y') | strcmp(tmp, 'yes') | isempty(tmp)
	%        labelblanks(songs_fid);
	%    end
end

fclose('all');

%-------------------------------------------------------------------------------

function filetype = get_filetype

global filetype

disp('What is type of sound file? [w]')
disp(' b = binary from mac')
disp(' w = wavefile (i.e. cbs)')
disp(' d = dcpfile')
disp(' f = foosong/gogo')
disp(' o = observer file (last/song channel)')
disp(' o1r = observer file (second to last)')


filetype = 'null';
while strcmp(filetype,'null')
  temp=input(' ','s');
  if (strcmp(temp,'b') | strcmp(temp,'B') | strcmp(temp,'bin') | strcmp (temp, 'binary'));
     filetype = 'b';
  elseif strcmp(temp,'w')  | strcmp(temp,'wav') | strcmp(temp,'W') | isempty(temp)
     filetype = 'w';
  elseif strcmp(temp,'d')  | strcmp(temp,'dcp') | strcmp(temp,'D')
     filetype = 'd';  
  elseif strcmp(temp,'f')  | strcmp(temp,'F')
     filetype = 'f';  
  elseif strcmp(temp,'o')  | strcmp(temp,'O') | strcmp(temp,'obs0r')
     filetype = 'obs0r';
  elseif strcmp(temp,'o1r')  | strcmp(temp,'O1r') | strcmp(temp,'obs1r')
     filetype = 'obs1r';  
  else
     disp('Unacceptable! Pick again')
  end
end

%-------------------------------------------------------------------------------

function [training_fid] = get_training_batchfile
% Prompts user for training batchfile

global train_path

training_fid = -1;
trainingfile = 0;
   disp('Select Training Batchfile');
   [trainingfile, train_path]=uigetfile('*','Select Training Batchfile');
   training_fid=fopen([train_path, trainingfile]);
   if training_fid == -1 | trainingfile == 0
        error(['Cannot open file: ' trainingfile]);
   end

%--------------------------------------------------------------------------
   
function [songs_fid] = get_labelling_batchfile
% Prompts user for batchfile of songs to be labelled

global label_path

songs_fid = -1;
songsfile = 0;
   disp('Select Batchfile of Songs to be Labeled');
   [songsfile, label_path]=uigetfile('*','Select Batchfile of Songs to be Labeled');
   songs_fid=fopen([label_path, songsfile]);
   if songs_fid == -1 | songsfile == 0
        error(['Cannot open file: ' songsfile]);
   end
   
%-------------------------------------------------------------------------------

function [sylls, index, labelIndex, id] = read_training(training_fid)
% Gets syllables and labels from training batchfile

global train_path
global filetype

sylls = [];
index =[];
labelIndex = [];
id=[];

while 1
   %get soundfile name
     soundfile = fscanf(training_fid,'%s',1);
   %end when there are no more notefiles 
     if isempty(soundfile)
        break
     end
     
   if (~(soundfile(1)=='/'))&(~strcmp(train_path,pwd))
       soundfile=[train_path,soundfile];
   end
   if (soundfile(end-4:end)=='.filt')   % if the batch file is of .filt files...
       notefile=[soundfile(1:end-5),'.not.mat'];
       filtfile=soundfile;
       if ~strcmp(filetype,'filt')
           disp('Filetype .filt found.')
       end
       filetype='filt';
   else
       if (soundfile(end-4:end)=='.cbin')  % if the batch file is of .cbin files...
           notefile=[soundfile, '.not.mat'];
           filtfile=[soundfile, '.filt'];
           if ~(strcmp(filetype,'obs0r')|strcmp(filetype,'obs1r'))
               disp('Filetype .cbin found.')
           end
           filetype='obs0r';
       else                              % if the files end in neither .cbin nor .filt
           notefile=[soundfile, '.not.mat'];
           filtfile=[soundfile, '.filt'];
           if ~(strcmp(filetype,'w'))
               disp('Assuming wave file format (file does not end in .cbin or .filt)')
           end
           filetype='w';
       end
   end

   % Skip file if it doesn't exist or has no .not.mat file
    if ~((exist(soundfile)|exist(filtfile)) & exist(notefile))   
        disp(['Skipping ', soundfile, ' (file or .not.mat file does not exist)'])
        continue
    end
    
    % Read filt file if it exists, otherwise read and filter raw song, saving a copy
    if exist(filtfile)
       [filtsong, fs]=read_filt(filtfile);
    else
        disp(['   Filtering ',soundfile])
        [song,fs]=soundin('', soundfile, filetype);
        filtsong=lt_bandpass(song,fs,300,8000);
        write_filt(filtfile, filtsong, fs);
    end
   
    load(notefile);
        
    %make sure labels is a column vector
    [x,y]=size(labels);
    if y > x
        labels=labels';
    end 
    
    disp(['Processing ' soundfile ' . . .']); 
    
    % convert onsets and offsets from ms to samples
    on = fix(onsets.*fs/1000);
    off = fix(offsets.*fs/1000);
    
    % normalize the amplitude of the syllables
%%    filtsong=normalize_sylls(filtsong, on, off);
    
    % For each syllable in the song file
    for i = 1:size(labels,1)
        
        % Skip certain labels
        if labels(i)=='0'|labels(i)=='-'
            continue
        end
           
        % Add syll
        sylls = [sylls; {filtsong(on(i):off(i))}];
        % Add label to index if it isn't there already
        if isempty(labelIndex) | isempty(find(labelIndex == labels(i)))
            labelIndex = [labelIndex; labels(i)];
        end
        % add label to index of sylls, in numerical form
        index = [index; find(labelIndex == labels(i))];
        id = [id; {soundfile},{onsets(i)}];
        
    end 
end

% Count number of examples and report results
counts = zeros(size(labelIndex));
for i = 1:size(labelIndex, 1)
    counts(i) = size(find(index == i),1);
    disp(['Found ' int2str(counts(i)) ' examples of ' labelIndex(i)]);
end
% Remove label if it has less than 8 examples
i = 1; inc=0;
while i <= size(labelIndex,1)
    if counts(i+inc) < 8
        disp(['Removing ' labelIndex(i)]);
        labelIndex = labelIndex(find(labelIndex ~= labelIndex(i)));
        sylls = sylls(find(index ~= i));
        index = index(find(index ~= i));
        for k = i:size(labelIndex,1)+1
            index(find(index == k)) = k - 1;
        end
        i = i - 1;
        inc=inc+1;
    else
       if counts(i+inc)<30
           disp(['Warning: ',labelIndex(i),' has only ',num2str(counts(i+inc)),' examples.']);
       end
    end
    i = i + 1;
end

%-------------------------------------------------------------------------------

function templates = make_templates(templateSylls)

% Calculates feature vectors for all example syllables

disp('Making templates.');

templates = [];

for i = 1:size(templateSylls,1)
    templates = [templates; feature_vect(templateSylls{i},32000)];
    if mod(i,100)==0                                                  
        disp([num2str(i),' of ',num2str(size(templateSylls,1)),' completed.'])
    end
end

%-------------------------------------------------------------------------------

function [template_cov,template_mean] = get_covariances(templates,index)

% Calculates the means and covariance matrices for each cluster of feature vectors

disp('Calculating covariance matrices.');

template_cov=[];
template_mean=[];
%
%  Commented out section is perfectly good without cyclic variables (e.g. with range 0-100, where 99-1 is -2, rather than 98)
%  New code intended to f
%
%for ii = 1:max(index)
%     template_cov=[template_cov;{cov(templates(index==ii,:))}];
%     template_mean=[template_mean;mean(templates(index==ii,:))];
%end

for ii = 1:max(index)
     these=templates(index==ii,:);
     meantemplate=mean(these);
     meantemplate(6)=50+100*angle(mean(exp(2i*pi*(these(:,6)-50)/100)))/(2*pi);   % Variables 6&7 are Cyclic
     meantemplate(7)=50+100*angle(mean(exp(2i*pi*(these(:,7)-50)/100)))/(2*pi);     
     template_mean=[template_mean;meantemplate];
     devs=these-meantemplate(ones(size(these,1),1),:);
     devs(:,6)=devs(:,6)+100*(devs(:,6)<-50)-100*(devs(:,6)>50);
     devs(:,7)=devs(:,7)+100*(devs(:,7)<-50)-100*(devs(:,7)>50);
     for l=1:length(meantemplate)
          for j=1:length(meantemplate)
               covtemplate(l,j)=0;
               for k=1:size(devs,1);
     	            covtemplate(l,j)=covtemplate(l,j)+devs(k,l)*devs(k,j);
               end
	       covtemplate(l,j)=covtemplate(l,j)/(size(devs,1)-1);
          end
     end
     template_cov=[template_cov;{covtemplate}];
end



%---------------------------------------------------------------------------------


function countblanks(songs_fid)

% Gives some concluding remarks. 

labeledCount = 0;
blankCount = 0;
numSongs = 0;

frewind(songs_fid);

while 1
   %get soundfile name
     songfile = fscanf(songs_fid,'%s',1);
   %end when there are no more notefiles 
     if isempty(songfile)
        break
     end

    numSongs = numSongs + 1; 
    notefile = [songfile '.not.mat'];
    if exist(notefile)
        load(notefile);
    else
        return
    end

    labels=labels';
    
    for i = 1:size(labels,1)
        if labels(i) == '-'
            blankCount = blankCount + 1;
        else
            labeledCount = labeledCount + 1;
        end
    end
end

disp(['Labeled ' int2str(numSongs) ' songs.']);
disp(['Gave labels to ' int2str(labeledCount) ' syllables.']);
disp(['Left ' int2str(blankCount) ' blank.']);

percentLabeled = fix(100*(labeledCount / (labeledCount + blankCount)));

disp([int2str(percentLabeled) ' percent of syllables labeled.']);

%--------------------------------------------------------------------------

function evaluate(templates,template_cov,template_mean,index)

% Evaluate the templates -- how tight are the syllable definitions, and how
% distinct are the different syllables.

global labelIndex

maxi=max(index);

disp('Evaluating templates...')
errors=zeros(maxi,maxi);
correct=zeros(maxi,1);

for k = 1:maxi
    kmean=template_mean(k,:);
    kdiff=(templates-kmean(ones(size(templates,1),1),:));
    kdiff(:,6)=kdiff(:,6)+100*(kdiff(:,6)<-50)-100*(kdiff(:,6)>50);
    kdiff(:,7)=kdiff(:,7)+100*(kdiff(:,7)<-50)-100*(kdiff(:,7)>50);
    kinvcov=((template_cov{k})^-1);
    d(:,k) = sqrt(sum((kdiff*kinvcov).*kdiff,2));
    disp([num2str(k),' of ',num2str(maxi),' clusters processed.']);
end
[tmp,class] = min(d');

for i = 1:length(templates(:,1))
%    d=zeros(maxi,1);
%    temp=templates(i,:);
%    for k = 1:maxi
%        templatek = templates(find(index == k),:);
%        d(k) = mahal(templates(i,:),templatek);
%    end
%    [tmp,class] = min(d');
    if class(i)==index(i)
        correct(class(i))=correct(class(i))+1;
    else
        errors(index(i),class(i))=errors(index(i),class(i))+1;
    end
    if mod(i,100)==0
        disp([num2str(i),' of ',num2str(length(templates(:,1))),' points processed.'])
    end
end

%for k = 1:maxi
%    kmean=template_mean(k,:);
%    kdiff=(featvect-kmean(ones(size(featvect,1),1),:));
%    kinvcov=((template_cov{k})^-1);
%    d(:,k) = sqrt(sum((kdiff*kinvcov).*kdiff,2));
%end


nf=length(templates(1,:));
distmat=zeros(maxi,maxi);
for i=1:maxi
    templatei = templates(find(index == i),:);
        kmean=template_mean(i,:); 
        kdiff=(templatei-kmean(ones(size(templatei,1),1),:));
        kdiff(:,6)=kdiff(:,6)+100*(kdiff(:,6)<-50)-100*(kdiff(:,6)>50);
        kdiff(:,7)=kdiff(:,7)+100*(kdiff(:,7)<-50)-100*(kdiff(:,7)>50);
        kinvcov=((template_cov{i})^-1);
        selfdist=sum(abs((kdiff*kinvcov).*kdiff),2);
        distmat(i,i)=mean(selfdist);
        disp([labelIndex(i),':   self-distance = ',num2str(mean(selfdist),3),' +/- ',num2str(std(selfdist),3)])
    for j=1:maxi
        if ~(j==i)
            templatej = templates(find(index == j),:);
            kdiff=(templatej-kmean(ones(size(templatej,1),1),:));
            kdiff(:,6)=kdiff(:,6)+100*(kdiff(:,6)<-50)-100*(kdiff(:,6)>50);
            kdiff(:,7)=kdiff(:,7)+100*(kdiff(:,7)<-50)-100*(kdiff(:,7)>50);
            dist=sum(abs((kdiff*kinvcov).*kdiff),2);
            dmean=mean(dist);
            distmat(i,j)=dmean;
            if (dmean<10*nf)
                disp(['     distance from ',labelIndex(j),' is only ',num2str(dmean,3),' +/- ',num2str(std(dist),3)])
            end
            if (errors(i,j)>0)
                disp(['    mislabeled as ',labelIndex(j),' in ',num2str(errors(i,j)),' of ',num2str(correct(i)+sum(errors(i,:))),' cases.'])
            end
        end
    end
    if ~any(errors(i,:))
        disp('    never mislabeled.')
    end
    disp(' ')
end

distmeds=median(distmat');

for i=1:maxi
    if (distmeds(i)<10)
        disp(['!!!!!!   Note ',labelIndex(i),' is very indistinct -- are you sure it is a single note?   (',num2str(distmeds(i),3),')']);
    end
    for j=1:maxi
        if (~(i==j))&(distmat(i,j)<5*nf)&(distmat(i,j)<distmeds(i)/2)
            disp(['!!!!!!!!   Are you certain that ',labelIndex(i),' and ',labelIndex(j),' are different syllables???'])
            disp(['   The distances between them are only ',num2str(distmat(i,j),3),' and ',num2str(distmat(j,i),3),'.'])
        end
    end
end

errnum=sum(sum(errors));
corrnum=sum(correct);
tot=errnum+corrnum;
disp(' ')
disp(['Net results: ',num2str(corrnum),' of ',num2str(tot),' labelled correctly.  (',num2str(errnum),' wrong).'])


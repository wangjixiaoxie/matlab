function lt_neural_v2_Finalize(DesiredClust, electrode_depth, Notes, LearningParams, deleteAudioFromMetadat)

% e.g. lt_neural_v2_Finalize(1, 1975, {'SUnit_0', 'Location_LMAN'}, {'29Apr2017-1200',{}})

%% LT 3/7/17 - by default removes the concatenated audio file from metadat 
% file size too large. can extract relatively easily from raw files.
% not done yet !!!! alternatively convert from double to single - confirmed
% that does not affect pitch estimates singifincatly (i.e. ADC is only
% 16bit)

%% LT 2/24/17 - Run in waveclus folder after happy with clustering -
% save information about data, metadata, etc. to a final structure that
% assists with final analysis.

% FOLDER NEEDS TO BE LIKE:
% /bluejay5/lucas/birds/bu77wh13/NEURAL/020317_LMANneural1/Chan9amp-BatchNight2150

%% === Inputs

% DesiredClust = []; % single number, which clust to save?
% electrode_depth = [];
% Notes = {};
    % HOW TO FORMAT NOTES
    % Notes = {'SUnit_0', 'Location_LMAN', 'Random stuff....'}; not single unit. in LMAN
    %   Location can be: 'LMAN', 'X', 'LMAN?', 'X?', '?'
    %   SUnit can be 0 or 1
    % e.g. lt_neural_v2_Finalize(1, 1900, {'SUnit_0', 'could get more song'})

% LearningParams{1} = '05Feb2017-1718'; % time of WN on [or when WN
% switched first time for this neuron]
% LearningParams{2} = {'07Feb2017-1509', '07Feb2017-1725'}; % cell aray,
% each entry one important date (e.g. wn direction change) - mainly for
% plotting


if ~exist('deleteAudioFromMetadat', 'var')
    deleteAudioFromMetadat = 0; % this removes field: metaDat.songDat, since is too large
end

%% save backup (overwrites old backup)

currdir = pwd;
targdir = ['/bluejay5/lucas/analyses/neural/'];
cd(targdir);
eval('!cp SummaryStruct.mat BACKUP/')
cd(currdir);


%% === query user as to which cluster to keep?

load('times_data')

max_clust_num = max(cluster_class(:,1));

if max_clust_num > 1
    assert(exist('DesiredClust', 'var')==1, 'PROBLEM - tell me what cluster num desired')
    assert(~isempty(DesiredClust), 'PROBLEM - tell me what cluster num desired')   
else
    DesiredClust =1;
end

assert(DesiredClust <= max_clust_num, 'PROBLEM - desired clust is too large')

  


%% get params

currdir = pwd;
slashes = findstr(currdir, '/');
dashes = findstr(currdir, '-');
uscores = findstr(currdir, '_');


birdname = currdir(slashes(4)+1:slashes(5)-1);
date_mmddyy = currdir(slashes(6)+1:slashes(6)+6);

tmp = find(uscores>slashes(6)+8 & uscores < slashes(7));
if isempty(tmp)
ExptID = currdir(slashes(6)+8:slashes(7)-1);
else
    ExptID = currdir(slashes(6)+8:uscores(tmp)-1);
end
batchfilename = currdir(dashes(end)+1:end);

Channel = currdir(slashes(end)+1:dashes(end)-1);
Channel = str2num(Channel(5:end-3));

disp(['birdname = ' birdname]);
disp(['date = ' date_mmddyy]);
disp(['ExptID = ' ExptID]);
disp(['batchfilename = ' batchfilename]);
disp(['Channel = ' num2str(Channel)]);
disp(['DesiredClust = ' num2str(DesiredClust)]);


%% Update metadat structure
disp('  ---  ');

clear SummaryStruct
load([targdir 'SummaryStruct']);


% -- does this bird exist?
birdindlogic = strcmp({SummaryStruct.birds.birdname}, birdname);
if ~any(birdindlogic) 
    % then bird doesn't exist
    SummaryStruct.birds(length(birdindlogic)+1).birdname = birdname;
    SummaryStruct.birds(length(birdindlogic)+1).neurons = [];
end

BirdInd = find(strcmp({SummaryStruct.birds.birdname}, birdname));

% -- how many neurons already exist for this bird?
if ~isfield(SummaryStruct.birds, 'neurons')
    SummaryStruct.birds.neurons = [];
end
numneurons = length(SummaryStruct.birds(BirdInd).neurons);
if numneurons ==0
    % if no neurons, then start from 1
    NeuronInd = 1;
else
    % if contains neurons, as whether this neuron already is entered.
    % if already entered, then replace data
    NeuronInd = find(strcmp({SummaryStruct.birds(BirdInd).neurons.dirname}, currdir) ...
        & [SummaryStruct.birds(BirdInd).neurons.clustnum]==DesiredClust); % match based on dir name + cluster number
    if isempty(NeuronInd)
        disp('New neuron! adding')
        % then add a new neuron
       NeuronInd = numneurons+1;
    else
        disp(['Old neuron (num ' num2str(NeuronInd) '), replacing old data']);
    end
end
   

%% ============ EXTRACT SOME USEFUL INFO

fid = fopen(batchfilename);
fline = fgetl(fid);

DatenumAll = [];
DatestrAll = {};
while ischar(fline)
   % time of start of each song file
[dtnum, datestring]=lt_neural_fn2datenum(fline);

DatenumAll = [DatenumAll dtnum];
DatestrAll = [DatestrAll datestring];

fline = fgetl(fid);
end


% ======= IS SINGLE UNIT?
inds = regexp(Notes, 'SUnit');

sunitstr = Notes{~cellfun(@isempty, inds)};
IsSU = str2num(sunitstr(inds{~cellfun(@isempty, inds)}+6));
if IsSU==1
    IsSU = 'yes';
else
    IsSU = 'no';
end


% =========== LOCATION.
inds = regexp(Notes, 'Location');
locationstr = Notes{~cellfun(@isempty, inds)};
probeLocation = locationstr(10:end);


%% ========== ENTER DATA
% KEEP EACH ISOLATED UNIT(CLUSTER) AS A SEPARATE NEURON

SummaryStruct.birds(BirdInd).neurons(NeuronInd).dirname = currdir;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).exptID = ExptID;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).batchfilename = batchfilename;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).date = date_mmddyy;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).channel = Channel;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).clustnum = DesiredClust;

SummaryStruct.birds(BirdInd).neurons(NeuronInd).electrode_depth = electrode_depth;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).Notes = Notes;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).NOTE_is_single_unit = IsSU;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).NOTE_Location = probeLocation;

SummaryStruct.birds(BirdInd).neurons(NeuronInd).Filedatestr_unsorted = DatestrAll;
SummaryStruct.birds(BirdInd).neurons(NeuronInd).Filedatenum_unsorted = DatenumAll;

if exist('LearningParams', 'var')
   SummaryStruct.birds(BirdInd).neurons(NeuronInd).LEARN_WNonDatestr = LearningParams{1};
   SummaryStruct.birds(BirdInd).neurons(NeuronInd).LEARN_WNotherImportantDates = LearningParams{2};
else
     SummaryStruct.birds(BirdInd).neurons(NeuronInd).LEARN_WNonDatestr = '';
   SummaryStruct.birds(BirdInd).neurons(NeuronInd).LEARN_WNotherImportantDates = {''};
end
    


%% ======== SAVE

save([targdir 'SummaryStruct'], 'SummaryStruct');
disp('DONE!! saved...')


%% ======== REMOVE AUDIO FROM METADAT (REDUCE FILE SIZE)

if deleteAudioFromMetadat==1
disp('Deleting songDat from MetaDat');    
    load('MetaDat.mat')
    
    metaDat = rmfield(metaDat, 'songDat');
    
    save(['MetaDat.mat'], 'metaDat', '-v7.3'); % save metadat
end


%% ======= MAKE NOTE IN DIRECTORY SAYING IS FINALIZED

tstamp = lt_get_timestamp(0);
fname = ['FINALIZED_' tstamp];

fid = fopen(fname, 'w');
fclose(fid);


%% ======= CONVERT FROM NEURONDATABASE TO SUMMARY STRUCT
if (0)
    
end


%% ======= DELETE LARGE RAW DATAFILES IF DESIRED



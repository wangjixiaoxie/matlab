function results =load_uisonganal_results(batch_file_filt)

% function
% results = load_uisonganal_results(batch_file_filt)
%
% This function prompts the user for a batchfile of song names and
% then a directory where the uisonganal processed data for these
% songs are stored (the results_dir).  This allows the user to use the same batch
% file that he used to analyze the songs using uisonganal.  This
% function then loads the .not.mat files (if they exist) from the
% results_dir.  The filename, labels,
% min_dur, min_int, offsets, onsets, sm_win and threshold for each
% file are put into data structs and returned in the data cell
% RESULTS.
% RESULTS{i}, for i=1:number of files, contains a struct with the
% following fields.
%   'filename'
%    'onsets'
%    'offsets'
%    'min_dur'
%    'min_int'
%    'labels'
%    'sm_win'
%    'threshold'
%
%  
% Example:  To access the name of the file for the 5th file loaded,
% type results{5}.filename   To access the onsets for this same
% file type results{5}.onsets. etc....
%
% The optional argument BATCH_FILE_FILT allows the user to input a
% file filter to speed the process of selection of the batch file.
% For instance if the batch file is called birdsong.batch and is in
% directory /net/bigmachine/mydir/birdsong/ one can set
% BATCH_FILE_FILT = '/net/bigmachine/mydir/birdsong/*.batch' and
% the gui will go to that directory for loading your batchfile.
%
% Dependancies:  filecb.m
%
% Note:  This code required that one run matlab with java turned ON.


% Get batchfile of song names.  This can be the same batch file you
% use for uisonganal.
if nargin==0
 [bfile, bpath]=uigetfile('*','Select batch file of song names');
else
  [bfile, bpath]=uigetfile(batch_file_filt,'Select batch file of songnames');
end

% Get initial directory.  
initdir = pwd;

% Get the results directory
btitle = sprintf('Dir to load results from batch file %s',bfile);
results_dir = filecb('browsedir',initdir,btitle);

bfid = fopen([bpath,bfile]);
ifile = 1;

while 1
    % To ** written by BD Wright.
    fname = fgetl(bfid);
    % Check for whitespace here!
    spaceflag = 0;
    if isspace(fname)
      spaceflag = 1;
    end
    if (~ischar(fname)) | isempty(fname) | spaceflag
      disp('End of batch file reached.');
      break
    end
    % **

    % get .not.mat filename 
    not_file = fullfile(results_dir,[fname,'.not.mat']);
    
    % Load this file (if it exists) and extract it's contents...
    if exist(not_file,'file')
      load(not_file);
      
      % fill data structures 
      
      tstruct.filename=fname;
      tstruct.onsets = onsets;
      tstruct.offsets = offsets;
      tstruct.min_dur = min_dur;
      tstruct.min_int = min_int;
      tstruct.labels = labels;
      tstruct.sm_win = sm_win;
      tstruct.threshold =threshold;
      
      % fill cell with struct
      results{ifile} = tstruct;
      ifile = ifile+1;
    else
      disp('no note file to load')
    end
  end

    
   len = length(results);
   
   lenstr = sprintf('The cell results contains %d structs. (One for each .not.mat file loaded.)',len);
   
 

 disp('################################################################################')
 disp(lenstr)
 disp('################################################################################')
 
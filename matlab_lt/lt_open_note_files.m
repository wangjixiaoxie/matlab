%% LT 2014 - use when in the bird directory
% arg = 0 --> opens 1_LABEL_NOTES and 1_NOTE_TO_SELF
% arg = 1 --> also tries to open 1) summary ppt, 2) analysis scripts

function lt_open_note_files(arg)

% open these 2 files stored in bird directory
try
    open 1_LABEL_NOTES.txt
catch err
    try
        open 1_LABELS_NOTES.txt
    catch err
        disp('starting new 1_LABEL_NOTES.txt');
        edit 1_LABEL_NOTES.txt
    end
end
try
    open 1_NOTE_TO_SELF.txt
catch err
    disp('starting new 1_NOTE_TO_SELF.txt');
    edit 1_NOTE_TO_SELF.txt
end

try
    open 1_ANALYSIS_LOG.txt
catch err
    disp('starting new 1_ANALYSIS_LOG.txt');
    edit 1_ANALYSIS_LOG.txt
end

try 
    open 1_LOG_EXPTS.txt
catch err
    disp('starting new 1_LOG_EXPTS.txt');
    edit 1_LOG_EXPTS.txt
end

try 
open 1_CONCLUSIONS.txt
catch err
    disp('starting new 1_CONCLUSIONS.txt');
    edit 1_CONCLUSIONS.txt
end
    




if(0)
    % try to open 1) summary ppt, and 2) analysis notes
    if arg==1;
        [birdname bluejaynum]=lt_get_birdname_date_from_dir(0);
        try
            string1_fn=['/home/lucast4/Desktop/NOTES/bird_summary/' birdname '_summary*']
            string1_fn=ls(string1_fn)
            string1=['!libreoffice -impress ' string1_fn];
            eval(string1);
        catch err
            disp('Summary ppt ERROR: likely does not exist, not in directory, or improperly named');
        end
        try
            string2 = ['/home/lucast4/MATLAB/analysis_lt/' birdname '_analysis_*'];
            string2_fn=ls(string2)
            open(string2_fn);
        catch err
            disp('Analysis string file ERROR: likely does not exist, not in directory, or improperly named');
        end
    end
end


% try to open analysis script
[birdname bluejaynum]=lt_get_birdname_date_from_dir(0);
analysis_files=dir(['/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/analysis_lt/' birdname '_analysis*']);

if length(analysis_files)>1; % then multiple analysis files
    for i=1:length(analysis_files);
        disp([num2str(i) '  -  ' analysis_files(i).name]);
    end
    
    kk=input('Type number of analysis file you wish to open (e.g. 1) ');
    
    edit(['/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/analysis_lt/' analysis_files(kk).name]);
    
elseif length(analysis_files)==1;
    edit(['/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/analysis_lt/' analysis_files.name]);
    
elseif length(analysis_files)==0;
    disp('No analysis script found');
end








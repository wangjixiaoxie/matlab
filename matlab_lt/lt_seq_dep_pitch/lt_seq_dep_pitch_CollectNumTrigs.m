function [TotalNumTrigs NotenumsUnique NoteNumsAllVect] = lt_seq_dep_pitch_CollectNumTrigs
%% for a given day goes through all songs and finds all triggers and counts them up.
% also outputs the notenum for those triggers (evtafv4)
lt_make_batch(3)
batchf = 'batch';

fid = fopen(batchf);
fline = fgetl(fid);

NumTrigsAll = [];
NoteNumsAll = {}; % one cell per song
NoteNumsAllVect = []; % all notes in one vector

disp('Counting trigs ...');

while ischar(fline)
    
%     disp(fline);
    
    % ============== CHECK RECF FOR TEMPLATE HITS
    recf = readrecf([fline(1:end-5) '.rec']);
   
    if isempty(recf)
    fline = fgetl(fid);         
        continue
    end
    
    numtrigs = length(recf.ttimes);
    
    assert(length(recf.pbname) == numtrigs, 'asdfasdf');
    NoteNums = [];
    for i=1:numtrigs
        tmp = strfind(recf.pbname{i}, 'Templ');
        notenum = recf.pbname{i}(tmp+8);
        NoteNums = [NoteNums notenum];
        NoteNumsAllVect = [NoteNumsAllVect notenum];
    end
    
    NumTrigsAll = [NumTrigsAll numtrigs];
    NoteNumsAll = [NoteNumsAll NoteNums];

    fline = fgetl(fid); 
end

TotalNumTrigs = sum(NumTrigsAll);
NotenumsUnique = unique(NoteNumsAllVect);

disp('DONE!');

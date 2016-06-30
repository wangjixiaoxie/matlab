%% LT 10/1 - added description
% goes through all files. outputs file and within file indices for when
% labeled syls occur
% INPUTS - should be self explanatory

% OUTPUTS
% vlsorfn - cell array of file names, one ind for each note detected
% vlsorind - position within file, one ind for each note  (vector)

function [vlsorfn, vlsorind] = jc_vlsorfn(batch,NOTE,PRENOTE,POSTNOTE)


ff=load_batchf(batch);
note_cnt = 0;

vlsorfn=cell(1);
vlsorind=[];
for ifn=1:length(ff)
    fn=ff(ifn).name;
    fnn=[fn,'.not.mat'];
    
    if (~exist(fnn,'file'))
        continue;
    end
    disp(fn);
    load(fnn);

    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE]);

    for ii = 1:length(p)
         note_cnt = note_cnt +1;       
         vlsorfn{note_cnt} = fn;
         vlsorind(note_cnt) = p(ii)+length(PRENOTE);
                
    end
     
end
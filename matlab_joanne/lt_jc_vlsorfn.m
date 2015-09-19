%% LT 11/1/14 - initiated outputs. so will not get error if the syl does not exist.
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
     
if isempty(vlsorind);
    disp('Syl not found in any song');
end
end
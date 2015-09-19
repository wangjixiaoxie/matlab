function [offsetlist] = modify_onsets(bt, NOTE, PRENOTE, POSTNOTE, templ, pretm, posttm, zeroval)
fs=32000;
%convert pretm and posttm to points. 

prepts=round(pretm*fs);
pstpts=round(posttm*fs);


ff=load_batchf(bt);
 offsetlist=[];
for ifn=1:length(ff)
   
    fn=ff(ifn).name;
    fnn=[fn];
    if (~exist(fnn,'file'))
        continue;
    end
    [data,fs]=ReadCbinFile(fnn);
    
    [sm,sp,t,f]=evsmooth(data,fs);

   
    %this code taken from findwnote4 is used to identify each occurence of
    %the note
    
    fnn=[fn,'.not.mat'];
    if (~exist(fnn,'file'))
        continue;
    end
    disp(fn);
    load(fnn);
    
    labels = lower(labels);
    labels(findstr(labels,'0'))='-';
   

    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE])+length(PRENOTE);
    mod_onset=onsets;    
    
    for ii = 1:length(p)
        if(length(onsets)==length(labels))
            ton=onsets(p(ii));
            %onsets are in ms, convert them to points.
            ton_pts=round((ton/1000)*32000)
            cur_data=sm((ton_pts-prepts):(ton_pts+pstpts));
            tmplen(1)=length(templ)
            tmplen(2)=length(cur_data)
            minlen=min(tmplen)
            templ=templ(1:minlen)
            cur_data=cur_data(1:minlen)
            %here call a function which will take as an input the template
            %and the data, and return an offset
            [offset]=calcoffset(templ, cur_data)
            offset=offset-zeroval;
            offsetlist=[offsetlist offset]
            %positive value of offset means time needs to be subracted.
            mod_onset(p(ii))=ton-(offset/fs)*1000;
        end
    end
    cmd=['save -append ' fnn ' mod_onset'];
    eval(cmd);
    
  
    end
    





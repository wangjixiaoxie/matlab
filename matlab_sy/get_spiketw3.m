%get_spiketw3 is being rewritten to deal for purposes of ra inactivation
%experiment
%TH=-2000;
bt='batch';
in_chans=2:3;
song_chan=1;
lpcnt=0;
wavesflag=0;
TH=get_thresh3(bt,song_chan,in_chans, in_chans)
filedata.params=[];

fid=fopen(bt,'r');
    
flag=0;
      
spktar=[];
wavesar=[];
spkampar=[];
timestamps=[];
        
while (1)
    fn=fgetl(fid);    
    if (~ischar(fn));
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    disp(fn);
    
    
    
    for jj=1:length(in_chans)
        spkt=[];
        spkamp=[];
        
        %if(~exist(teststr))
                [data,fs,spkind,spkamp,waves,tmstmps]=tetanaltw2(fn,TH(jj),song_chan,in_chans(jj),wavesflag);
                if(flag==0)
                    off_lng=ceil(length(data(:,1))/fs); 
                    flag=1;
                end
                
                spkt=spkind/fs;
               % spkt=spkt+off_lng*lpcnt;
              

                
                    clustspk(jj).spkt=spkt;
                    clustspk(jj).spkamp=spkamp; 
                
    end
    eval(['save ',fn,'.clust clustspk']);
       end
    
    fclose(fid);
    
 
 
 %every single file gets its own clustspk structure saved
 
 
 

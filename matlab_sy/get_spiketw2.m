%clear all;
%This was rewritten to work with mattias' clustering program.



%get_thresh will find an area two standard deviations above the noise level


%TH=-2000;
tet_chans=2:5;
song_chan=1;
lpcnt=0;
wavesflag=1;
bt='batch1';

TH=get_thresh(bt,song_chan,tet_chans)
if(wavesflag)
    paramnames={'spktime' 'spkamp1' 'spkamp2' 'spkamp3' 'spkamp4' 'spkwid1' 'spkwid2' 'spkwid3' 'spkwid4'};
else
    paramnames={'spktime' 'spkamp1' 'spkamp2' 'spkamp3' 'spkamp4'}
end
    site='site2'

    flag=0;
fid=fopen(bt,'r');

spktar=[];
wavesar=[];
spkampar=[];
spkwidar=[];
timestamps=[];
spkwid1=[];
spkwid2=[];
while (1)
	%each iteration through this loop is a single file
    fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	disp(fn);
    teststr=strcat(fn,',spk');
   %if(~exist(teststr))
        [data,fs,spkind,spkamp,waves,tmstmps]=tetanaltw2(fn,TH,song_chan,tet_chans,wavesflag);
        if(flag==0)
           off_lng=ceil(length(data(:,1))/fs); 
           flag=1;
        end
        waves=int16(waves);
        spkt=spkind/fs;
        spkt=spkt+off_lng*lpcnt;
        spktar=[spktar;spkt];
        spkampar=[spkampar;spkamp];
        clear spkwid
        clear spkwid2
        if(wavesflag)
           
           spkwid1=Find_Spk_Width2(waves);
               
            for i=1:4
                spkwid2(:,i)=spkwid1(:,1,i);
            end
            spkwidar=[spkwidar;spkwid2];
            wavesar=cat(1,wavesar,waves);
            timestamps=[timestamps;tmstmps];
        end
    lpcnt=lpcnt+1;

end
fclose(fid);
 if(wavesflag)
    filedata.params=[spktar spkampar spkwidar];
 else
    filedata.params=[spktar spkampar];
 end
 filedata.paramnames=paramnames;
 filedata.filename=site;
 filedata.fileloc=pwd;
 outstr=[site '_params'];
 eval(['save ',outstr,' filedata']);
 if (wavesflag)
    waves=wavesar;
    eval(['save ',site,'.mat timestamps waves']);
 end
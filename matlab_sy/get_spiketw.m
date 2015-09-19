%clear all;
sfct=(1e3/2^15);
tet_chans=2:5;
song_chan=1;
bt='batch';
%TH=get_thresh(bt,song_chan,tet_chans);
%TH=-TH/sfct
TH=-5000;
fid=fopen(bt,'r');
while (1)
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
        [data,fs,spkind,spkamp]=tetanal(fn,TH,song_chan,tet_chans);
        %waves=int16(waves);
        spkt=spkind/fs;
        %spkt=uint32(spkt);
        eval(['save ',fn,'.spk spkt spkamp']);
    %end
    end
fclose(fid);

%clear all;
sfct=(1e3/2^15);
in_chans=2:3;
song_chan=1;
bt='batch';
TH=get_thresh3(bt,song_chan,in_chans,in_chans);
%TH=-TH/sfct
%TH=-5000;
fid=fopen(bt,'r');
wavesflag=0;
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
    for ii=1:length(in_chans)
        inchn=in_chans(ii);
    %if(~exist(teststr))
        [data,fs,spkind,spkamp]=tetanaltw2(fn,TH(jj),song_chan,inchn,wavesflag);
        %waves=int16(waves);
        spkt=spkind/fs;
        %spkt=uint32(spkt);
        teststr=strcat(fn,num2str(inchn));
        eval(['save ',teststr,'.spk spkt spkamp']);
    %end
    end
    end
fclose(fid);

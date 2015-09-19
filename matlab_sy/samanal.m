%useful for analyzing sam's data
batch='/doyale1/twarren/pu3bu86/tw_data_vols/batch.keep.anal';
chn2plot=2;
song_chan=1'
tet_chans=2:6;
TH=[0 0];
fnm={};
fid=fopen(batch,'r');
upper=TH(2);
thr=TH;
samplerate=32000;
chmast=[];
indmast=[];
cnt=0;
while (1)
	fn=fgetl(fid)
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
        continue;
	end
	
	cnt=cnt+1;
	fnm{cnt}=fn;
end
fclose(fid);

for ifn=1:length(fnm);

    fn=fnm{ifn};
    %run tetanal on the file get spike times and spkamps
	
    rd=readrecf(fn);
    %read in data
    if(isempty(rd)==0)
        [data,fs]=Readcbinfile(fn);
    end
    
    sfct=(1e3/2^15);
    data=-sfct.*data;
    snd=data(:,1)';
    ch=data(:,chn2plot)';
    subplot(2,1,1)
    plot(snd);
    title(fn);
    subplot(2,1,2)
    plot(ch,'b');
    figure;
    
    [thr]=get_thresh2(data,song_chan,tet_chans,chn2plot,TH,fs,fn)
    
    chmast=[chmast ch];
    ind=ifn*ones(size(ch));
    %indmast=[indmast;ind];

    close all;
    [outstruct(ifn).spktimes,ind]=analbird(upper,thr,snd,ch,samplerate);
    outstruct(ifn).fn=fn;
    %outstruct(ifn).ttimes=ttimes;
end
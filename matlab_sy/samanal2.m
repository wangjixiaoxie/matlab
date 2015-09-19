%useful for analyzing sam's data
batch='batch';
%batch='batch'
chn2plot=3;
song_chan=1;
tet_chans=2:5;
TH=[0 0];
indcount=1;
fid=fopen(batch,'r');
upper=TH(2);
thr=TH(1);
samplerate=32000;
chmast=[];
indmast=[];
indval=[];
cnt=0;
fnm=[];
while (1)
	fn=fgetl(fid);
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
    outdata(ifn).fn=fn;
    %run tetanal on the file get spike times and spkamps
	
    rd=readrecf(fn);
    %read in data
    if(isempty(rd)==0)
        [data,fs]=Readcbinfile(fn);
    end
    
    outdata(ifn).tbefore=rd.tbefore;
    outdata(ifn).tafter=rd.tafter;
    outdata(ifn).ttimes=rd.ttimes;
    outdata(ifn).nsamp=rd.nsamp
    sfct=(1e3/2^15);
    data=-sfct.*data;
    snd=data(:,1)';
    ch=data(:,chn2plot)';
    
    [TH]=get_thresh2(data,song_chan,tet_chans,chn2plot,TH,fs,fn)
    
    chmast=[chmast ch];
    indval=[indval indcount];
    indcount=indcount+length(data);

end
    indval=[indval indcount];
close all;
    [ind]=analbird2(TH(2),TH(1),snd,chmast,samplerate);
    for(ifn=1:length(fnm))
        sel_ind=find(ind<indval(ifn+1)&ind>indval(ifn))
        outdata(ifn).spktms=(ind(sel_ind)-indval(ifn))/32000
        
    end
%NAME A BATCHFILE
%batchf='batch';
fid =fopen(batchf,'r');
ff=[];
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	ff(length(ff)+1).name=fn;
end
fclose(fid);

logfid=fopen('LOGFILE','a');

if (exist('BREAKDOWN')~=7)
	eval(['!mkdir BREAKDOWN']);
end

%set these params
NEWTAFTER=1; % new time of silence to calose a data file
song_chan=2; % which channel has the sound in it
TH=0.17;     %raw sound threshold
%%%%%%%%%%%%%%%%%%%% set above params
ii=0;
while (ii<length(ff))
	ii=ii+1;

	fn=ff(ii).name;
	disp(fn);
	rd=readrecf(fn);
	nchan=rd.nchan;
	nsamp=rd.nsamp;
	tbefore=rd.tbefore;
	tafter=rd.tafter;
	fs=rd.adfreq;

	%if (tafter<2.9)
	%	continue;
	%end
	%if (tbefore<2.9)
	%	continue;
	%end

	tmp=findstr(fn,'.');EXT=fn(tmp(end):end);
	filecnt=fn((tmp(end-1)+1):(tmp(end)-1));
	ll=length(filecnt);
	filecnt=str2num(filecnt)*100;

	tmpfilen=[fn(1:tmp(end)),'tmp'];
	tmpdata=load(tmpfilen).';
	nn=length(rd.thresh);
	tmpdata=reshape(tmpdata,[nn,length(tmpdata)./nn]).';
	tmpdatat=[0:size(tmpdata,1)-1]*128*2/fs+rd.tbefore;

	[data,fs]=evsoundin('',fn);


	% first pull out tbefore amount of data
	nbuf=ceil(fs*tbefore);
	buf=zeros([nbuf,nchan]);
	nafter=ceil(fs*NEWTAFTER);
	pp=find(abs(data(:,song_chan))>=TH);
	pp=[pp(find((pp>nbuf)&(pp<(size(data,1)-ceil(tafter*fs)))));length(data)];

	ppp=find(diff(pp)/fs >= NEWTAFTER);
	fileindx=nbuf+1;
	buf=data(1:nbuf,:);
	for jj=1:length(ppp)
		if (fileindx-nbuf<1)
			disp(['hey : ',fn,' ',num2str(jj),' fi=',num2str(fileindx)]);
		end
		inds=[max(fileindx-nbuf,1):min(pp(ppp(jj))+nafter,size(data,1))];
		filecnt=filecnt+1;
		tmp=findstr(fn,'.');
		newfn=['BREAKDOWN/',fn(1:tmp(end-1)),num2str(filecnt)];
		fid=fopen([newfn,EXT],'w','b');
		fwrite(fid,reshape(data(inds,:).',[numel(data(inds,:)),1]).','float');
		fclose(fid);
		rdtemp=rd;
		rdtemp.tafter=NEWTAFTER;
		rdtemp.nsamp=length(inds);
		tt=rdtemp.ttimes;
		tt=tt(find((tt>=(inds(1)+nbuf)*1e3/fs)&(tt<=(inds(end)-nafter)*1e3/fs)));
		rdtemp.ttimes=tt-(inds(1)*1e3/fs);
		wrtrecf([newfn,'.rec'],rdtemp);


		[y,stind]=min(abs(tmpdatat-((inds(1)+nbuf)/fs)));
		[y,enind]=min(abs(tmpdatat-((inds(end)-nafter)/fs)));
		tmp=[stind:enind];
		tmp2=tmpdata(tmp,:);
		tmp2=reshape(tmp2.',[numel(tmp2),1]);
		fid=fopen([newfn,'.tmp'],'w');
		for kk=1:length(tmp2)
			fprintf(fid,'%.5f ',tmp2(kk));
		end
		fclose(fid);

		if (exist([fn,'.not.mat']))
			load([fn,'.not.mat']);
			tmp=find((onsets>=(inds(1)*1e3/fs))&...
                                (offsets<=(inds(end)*1e3/fs)));
			labels=labels(tmp);
			onsets=onsets(tmp)-inds(1)*1e3/fs;
			offsets=offsets(tmp)-inds(1)*1e3/fs;
			eval(['save ',newfn,EXT,'.not.mat Fs labels min_dur min_int offsets onsets sm_win threshold']);
		end
			
		fprintf(logfid,'%s :\n',newfn);
		fprintf(logfid,'%d \n',fileindx);
		fprintf(logfid,'%.5e %.5e \n',inds(1)/fs,inds(end)/fs);
		fprintf(logfid,'%d - %d \n',stind,enind);
		fprintf(logfid,'\n');
		
		fileindx=pp(ppp(jj)+1);
	end
	tmp=findstr(fn,'.');
	%delete([fn(1:tmp(end)),'*']);
end
fclose(logfid);

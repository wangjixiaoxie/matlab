%NAME A BATCHFILE
%batchf='batch';
ff=load_batchf(batchf);

logfid=fopen('LOGFILE','a');

if (exist('BREAKDOWN')~=7)
	eval(['!mkdir BREAKDOWN']);
end

%set these params
TAFTER =1; % new time of silence to calose a data file
TBEFORE=2; % new time of silence to calose a data file
TH=0.005;  %raw sound threshold
%bird_name='bk12w65';
year='07';
%%%%%%%%%%%%%%%%%%%% set above params
ii=0;filecnt=0;
while (ii<length(ff))
	ii=ii+1;

	fn=ff(ii).name;
	disp(fn);

	[data,fs]=wavread(fn);

	data2=bandpass(data,fs,500,1e4,'hanningfir');
	data=data2;clear data2

	rd=readrecf(fn);

	tmp=findstr(fn,'_');
	mon=fn(tmp(2)+1:tmp(3)-1);
	if length(mon)<2; mon=['0',mon];end

	dy=fn(tmp(3)+1:tmp(4)-1);
	if length(dy)<2; dy=['0',dy];end

	hr=fn(tmp(4)+1:tmp(5)-1);
	if length(hr)<2; hr=['0',hr];end

	mn=fn(tmp(5)+1:tmp(6)-1);
	if length(mn)<2; mn=['0',mn];end

	tmp2=findstr(fn,'.wav');
	sc=fn(tmp(6)+1:tmp2(1)-1);
	if length(sc)<2; sc=['0',sc];end

	%get an evtaf style name
	fn_new_base=[bird_name,'_',dy,mon,year,'_',hr,mn,sc,'.'];

	% first pull out tbefore amount of data
	nbuf = ceil(fs*TBEFORE);
	buf  = zeros([nbuf,1]);
	nafter = ceil(fs*TAFTER);
	

	if (length(data)<=nbuf+nafter)
		filecnt=filecnt+1;
		wavwrite(data,fs,16,['BREAKDOWN/',fn_new_base,num2str(filecnt),'.wav']);
		continue;
	end

	pp=find(abs(data)>=TH);
	pp=[pp(find((pp>nbuf)&(pp<(length(data)-ceil(TAFTER*fs)))));length(data)];

	ppp=find(diff(pp)/fs >= TAFTER);
	fileindx=nbuf+1;
	buf=data(1:nbuf);
	for jj=1:length(ppp)
		if (fileindx-nbuf<1)
			disp(['hey : ',fn,' ',num2str(jj),' fi=',num2str(fileindx)]);
		end
		inds=[max(fileindx-nbuf,1):min(pp(ppp(jj))+nafter,length(data))];
		filecnt=filecnt+1;
		newfn=['BREAKDOWN/',fn_new_base,num2str(filecnt),'.wav'];
		wavwrite(data(inds),fs,16,newfn);

		%if (exist([fn,'.not.mat']))
		%	load([fn,'.not.mat']);
		%	tmp=find((onsets>=(inds(1)*1e3/fs))&...
                %               (offsets<=(inds(end)*1e3/fs)));
		%	labels=labels(tmp);
		%	onsets=onsets(tmp)-inds(1)*1e3/fs;
		%	offsets=offsets(tmp)-inds(1)*1e3/fs;
		%	eval(['save ',newfn,EXT,'.not.mat Fs labels min_dur min_int offsets onsets sm_win threshold']);
		%end
			
		fprintf(logfid,'%s :\n',newfn);
		fprintf(logfid,'%d \n',fileindx);
		fprintf(logfid,'%.5e %.5e \n',inds(1)/fs,inds(end)/fs);
		fprintf(logfid,'\n');
		
		fileindx=pp(ppp(jj)+1);
	end
	tmp=findstr(fn,'.');
	%delete([fn(1:tmp(end)),'*']);
end
fclose(logfid);

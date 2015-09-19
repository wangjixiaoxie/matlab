function mk_tempf(batchfile,templ,PREDATATIME,CHANSPEC,EVMODE);
%mk_tempf(batchfile,templ,PREDATATIME,CHANSPEC,EVMODE);

if (~exist('EVMODE'))
	EVMODE=0;
elseif (length(EVMODE)==0)
	EVMODE=0;
end
if ~exist('CHANSPEC')
    CHANSPEC='obs0';
end
if (~exist('PREDATATIME'))
    PREDATATIME=2.0;
end

NTEMPL=size(templ,2);
fid=fopen(batchfile,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));break;end;
    if (~exist(fn,'file'))
        continue;
    end
    disp(fn);
    [dat,fs]=evsoundin('',fn,CHANSPEC);

    dat=dat(fix(PREDATATIME*fs):end);
    vals=evtafsim(dat,fs,templ,EVMODE);

    tt=[1:length(vals)]*2*size(templ,1)/fs;

    %pp=find(tt>=PREDATATIME);
    %vals=vals(pp,:);

    vals2=zeros([numel(vals),1]);
    cnt=1;
    for ii=1:size(vals,1)
        for jj=1:size(vals,2)
            vals2(cnt)=vals(ii,jj);
            cnt=cnt+1;
        end
    end
    vals=vals2;
    clear vals2;
    
	ppp=findstr(fn,'.cbin');
	if (length(ppp)==0)
		ppp=findstr(fn,'.ebin');
	end
	fnt=[fn(1:ppp-1),'X.tmp'];
	fid2=fopen(fnt,'w');
	for jj=1:length(vals)
        	fprintf(fid2,'%.5e\n',vals(jj));
	end
	fclose(fid2);
end
fclose(fid);
return;

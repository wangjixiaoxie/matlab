function mk_tempf(batchfile,templ,PREDATATIME,CHANSPEC);
%mk_tempf(batchfile,templ,PREDATATIME,CHANSPEC);

if ~exist('CHANSPEC')
    CHANSPEC='obs0';
end
if (~exist('PREDATATIME'))
    PREDATATIME=4.0;
end

NTEMPL=size(templ,2);
fid=fopen(batchfile,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));break;end;
    if (~exist(fn,'file'))
        disp('hey');
        continue;
    end
    disp(fn);
    [dat,fs]=evsoundin('',fn,CHANSPEC);

    vals=evtafsim(dat,fs,templ);
    tt=[1:length(vals)]*256/fs;

	pp=find(tt>=PREDATATIME);
	vals=vals(pp);
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
	fnt=[fn(1:ppp),'tmp'];
	fid2=fopen(fnt,'w');
	for jj=1:length(vals)
        fprintf(fid2,'%.5e\n',vals(jj));
	end
	fclose(fid2);
end
fclose(fid);
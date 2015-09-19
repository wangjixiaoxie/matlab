function label_trigs(batch,NT,CSPEC,threshold,OVERWRT,USEX,min_int,min_dur,sm_win,filter_type);
%label_trigs(batch,NT,CSPEC,threshold,OVERWRT,USEX,min_int,min_dur,sm_win,filter_type);
% every trigger is labeled as NT
% if OVERWRT==0 then it will skip files with .not.mat already present
%
% label_trigs('batch.rand','a','obs0',threshold,0,0,min_int,min_dur);
% label_trigs('batch.rand','a','obs0',threshold,0,0,min_int,min_dur,2,'hanningfirff');
% label_trigs('batch.catch','a','obs0',threshold,0,0,min_int,min_dur,2,'hanningfirff');

%

if (~exist('filter_type','var'))
    filter_type='hanningfir';
end
if (~exist('USEX','var'))
    USEX=0;
end
if (~exist('OVERWRT','var'))
    OVERWRT=1;
end

if (~exist('sm_win','var'))
    sm_win=2.0;
elseif (length(sm_win)==0)
    sm_win=2.0;
end

if (~exist('min_int','var'))
    min_int=5.0;
end

if (~exist('min_dur','var'))
    min_dur=30.0;
end

if (exist('CSPEC','var'))
    if (~exist('threshold'))
        disp(['Need to specify seg thresh!']);
        return;
    end
end

fid=fopen(batch,'r');
if (fid==-1)
	disp(['could not find batch file :',batch]);
	return;
end

while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end

	if (~exist(fn,'file'));
		continue;
	end
    if (~exist([fn,'.not.mat'],'file'))
        [dat,Fs]=evsoundin('',fn,CSPEC);
        sm=SmoothData(dat,Fs,1,filter_type);
        [ons,offs]=SegmentNotesJC(sm,Fs,min_int,min_dur,threshold);
        onsets=ons*1e3;offsets=offs*1e3;
        fname=fn;
        labels=char(ones([1,length(onsets)])*fix('-'));
    else
        if (OVERWRT==0)
            continue;
        end
        load([fn,'.not.mat']);
    end
    
    disp(fn);

	rd=readrecf(fn,USEX);
	if (~isfield(rd,'ttimes'))
		rd.ttimes=[];
	end
	for ii = 1:length(rd.ttimes)
		pp=find((onsets<=rd.ttimes(ii))&(offsets>=rd.ttimes(ii)));
		if (length(pp)>0)
			labels(pp(end))=NT;
		end
	end

	fname=fn;
	cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
	                      'offsets onsets sm_win threshold'];
	eval(cmd);
    clear sm;
end
fclose(fid);
return;

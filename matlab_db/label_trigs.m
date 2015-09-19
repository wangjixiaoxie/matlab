function label_trigs(batch,NT,CSPEC,threshold,OVERWRT,USEX,min_int,min_dur,sm_win);
%label_trigs(batch,NT,CSPEC,threshold,OVERWRT,USEX,min_int,min_dur,sm_win);
% every trigger is labeled as NT
% if OVERWRT==0 then it will skip files with .not.mat already present
%

if (~exist('USEX'))
    USEX=0;
end
if (~exist('OVERWRT'))
    OVERWRT=1;
end

if (~exist('sm_win'))
    sm_win=2.0;
end

if (~exist('min_int'))
    min_int=5.0;
end

if (~exist('min_dur'))
    min_dur=30.0;
end

if (exist('CSPEC'))
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
    if(isempty(readrecf(fn)))
        continue
    end

    if (~exist([fn,'.not.mat'],'file'))

        [dat,Fs]=evsoundin('',fn,CSPEC);
        if(size(dat)==0)
            continue
        end
        sm=evsmooth(dat,Fs,0);
        [ons,offs]=SegmentNotes(sm,Fs,min_int,min_dur,threshold);
        onsets=ons*1e3;offsets=offs*1e3;
        fname=fn;
        labels=char(ones([1,length(onsets)])*45);
    else
        if (OVERWRT==0)
            continue;
        end
        load([fn,'.not.mat']);
    end

    %disp(fn);

        rd=readrecf(fn,USEX);
        if (~isfield(rd,'ttimes'))
                rd.ttimes=[];
        end
        for ii = 1:length(rd.ttimes)
                pp=find((onsets<=rd.ttimes(ii))&(offsets>=rd.ttimes(ii)));
                if (length(pp)>0)
                        labels(pp(end))=NT;
                end
       

        fname=fn;
        cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
                              'offsets onsets sm_win threshold'];
        eval(cmd);
    clear sm;
    end
end
fclose(fid);
return;
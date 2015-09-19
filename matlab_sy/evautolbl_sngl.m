function outdat=evautolbl_sngl(batch,TEMPLNT,SPTH,NTTH,NOTE,DOLABEL,CHANSPEC,THRESH,OVERLBL);
% outdat=evautolabel(batch,TEMPLNT,SPTH,NTTH,NOTE,DOLABEL,CHANSPEC,THRESH,OVERLBL);
%templnt is the note template
%SPTH is the power level in the FFT set which is set to -1 below 1 above
%NTTH is the level needed for a positive identification
%NOTE is the character to albel it
%DOLABEL =1 actually write it out
%OVERLBL = 0 will not label song files that have .not.mat present
%        = 1 will label everything

if (~exist('DOLABEL'))
    DOLABEL = 0;
end

if (~exist('OVERLBL'))
    OVERLBL=0;
end

if (~exist('CHANSPEC'))
    CHANSPEC='obs0r';
    disp(['ChanSpec = ',CHANSPEC]);
end

cnt=0;
p1=find(TEMPLNT>=SPTH); p2 = find(TEMPLNT<SPTH);
TEMPLNT(p1) = 1;TEMPLNT(p2) = -1;
fid=fopen(batch,'r');
fid2=fopen([batch,'.evlabeled'],'w');
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end

    fnnm=[fn,'.not.mat'];
    [tmpn2,tmpn,ext]=fileparts(fn);
    if (exist(fnnm,'file'))
        if (OVERLBL==0)&(DOLABEL==1)
            continue;
        end
        load(fnnm);
        if (strcmp(ext,'.ebin'))
            [dat,Fs]=readevtaf(fn,'0r');
        else
            [dat,Fs]=evsoundin('',fn,CHANSPEC);
        end
        [sm,sp,t,f]=evsmooth(dat,Fs,0.0);sp=abs(sp);
    else
        if (strcmp(ext,'.ebin'))
            [dat,Fs]=readevtaf(fn,'0r');
        else
            [dat,Fs]=evsoundin('',fn,CHANSPEC);
        end
        [sm,sp,t,f]=evsmooth(dat,Fs,0.0);sp=abs(sp);
        if (exist('THRESH'))
            threshold = THRESH;
        else
            threshold = 5e3;
        end
        min_int = 5;
        min_dur = 30;
        sm_win=2;
        [onsets,offsets]=SegmentNotes(sm,Fs,min_int,min_dur,threshold);
        onsets = onsets*1e3;offsets=offsets*1e3;
        labels = char(ones([1,length(onsets)])*48);
    end
    disp(fn);

    nt = zeros(size(TEMPLNT));
    tlen = size(TEMPLNT,2);
    dt = t(2)-t(1);
    tmpvec=zeros([length(onsets),1]);
    for ii = 1:length(onsets)
        onind = floor(onsets(ii)*1e-3/dt) + 1;
        ofind = floor(offsets(ii)*1e-3/dt) + 1;
        ntlen = (ofind-onind+1);
        if (ntlen > size(TEMPLNT,2))
            ofind = onind + size(TEMPLNT,2) - 1;
            ntlen = (ofind-onind+1);
        else
            % REMOVE
            ntlen = size(TEMPLNT,2);
            ofind = onind + ntlen-1;
        end
        nt = 0.0*nt;
        if ((ofind<size(sp,2))&(onind<size(sp,2)))
            nt(:,1:ntlen) = sp(1:size(TEMPLNT,1),onind:ofind);
            %nt = sp(1:size(TEMPLNT,1),onind:ofind);
            p1 = find(nt>=SPTH); p2 = find(nt<SPTH);
            nt(p1) = 1; nt(p2) = -1;
            match = sum(sum(nt.*TEMPLNT));
            tmpvec(ii) = match;
        end
    end
    cnt=cnt+1;
    outdat(cnt).fn=fn;
    outdat(cnt).mvals = tmpvec;
    outdat(cnt).onsets= onsets;
    outdat(cnt).offsets = offsets;
    outdat(cnt).labels = labels;
    pp = find(tmpvec>=NTTH);
    %labels = char(ones([1,length(onsets)])*48);
    for ii = 1:length(pp)
        labels(pp(ii))=NOTE;
    end
    if (DOLABEL==1)
        cmd = ['save ',fnnm,...
            ' Fs labels min_dur min_int offsets onsets sm_win threshold'];
        eval(cmd);
        fprintf(fid2,'%s\n',fn);
    end
end
fclose(fid);
fclose(fid2);
return;

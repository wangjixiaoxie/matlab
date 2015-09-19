function [vals,trigs]=jc_triglabel(batch,NOTE,UseDotTrigFile,DoLowerCase,NoDoCatch,ADDX);
% [vals,trigs]=triglabel(batch,NOTE,UseDotTrigFile,DoLowerCase,NoDoCatch,ADDX);
% normal batch file
% USEDOTTRIGFILE = 1 means to use the .rec file trigger times
%vals = [Nhit,Nnt,Ntrig]
% trigs.fn=fn;
% trigs.ftime=hr;
% trigs.ttimes=tt;
% trigs.ftime   = time of the file
% trigs.toffset = trigger offset from note onset
% trigs.tnotes  = labels where triggers happened
% trigs.totalnt = total # of notes
% trigs.totaltargetnt = total # of target notes
% trigs.nmatch = # of hits
% trigs.nontrignt  = non target triggered notes
% trigs.trigmintmp = evtaf trigger valuse min during the note
% trigs.ntmintmp = ntmintmp
%

if (~exist('NoDoCatch'))
	NoDoCatch=0;
end

if (~exist('ADDX'))
	ADDX=0;
end

if (~exist('UseDotTrigFile'))
    UseDotTrigFile=0;
else
	if (length(UseDotTrigFile)<1)
		UseDotTrigFile=0;
	end
end

if (~exist('DoLowerCase'))
    DoLowerCase=0;
end

fid=fopen(batch,'r');
cnt=0;
while (1)
    fn = fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    fnm=[fn,'.not.mat'];
    if (~exist(fnm,'file'))
        disp(['NOTMAT does not exist']);
        continue;
    end
    disp(fn);
    load(fnm);
    if (length(findstr(labels,'0'))==length(labels))
	    disp(['Hey']);
	    continue;
    end

    if (NoDoCatch==1)
	    rdtmp=[];
	    rdtmp=readrecf(fn);
	    if (length(rdtmp)>0)
		    if (rdtmp.iscatch)
			    continue;
		    end
	    end
    end
    
    if (DoLowerCase==1)
        pos = findstr(lower(labels),lower(NOTE));
    else
        pos = findstr(labels,NOTE);
    end
    [pth,tnm,ext]=fileparts(fn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn,'0r');
    else
        [dat,fs]=evsoundin('',fn,'obs0r');
    end
    if (UseDotTrigFile==0)
        if (strcmp(ext,'.ebin'))
            [trg,fs]=readevtaf(fn,'1r');
            bnd=0.5;
        else
            [trg,fs]=evsoundin('',fn,'obs1r');
            bnd=1e4;
        end
        trg=trg-min(trg);trg=trg./max(trg);bnd=0.5;
        tt=find((trg(1:end-1)<=bnd)&(trg(2:end)>bnd));
        tt=tt/fs;
    else
        rd=readrecf(fn);
	if (~isfield(rd,'ttimes'))
		rd.ttimes=[];
	end
        tt = rd.ttimes/1e3;
    end
    fpp = findstr(fn,ext);
    if (ADDX==1)
    	fnt = [fn(1:fpp(end)-1),'X.tmp'];
    else
    	fnt = [fn(1:fpp(end)),'tmp'];
    end
    if (exist(fnt,'file'))
        tdata=load(fnt);
        tdatat = [0:length(tdata)-1]*256/fs + rd.tbefore;
    else
	    tdatat=[];
    end
    cnt=cnt+1;
    repvals = get_repeats(labels);

    hr=fn2date(fn);
    trigs(cnt).fn=fn;
    trigs(cnt).ftime=hr;
    trigs(cnt).ttimes=tt;
    nt=[];nfndcnt=0;mintmp=[];toffset=[];ttilend=[];repcnt=[];
    trigind=[];
    for ii = 1:length(tt)
        ptmp=find((onsets<=tt(ii)*1e3)&(offsets>=tt(ii)*1e3));
        if (length(ptmp)>0)
            if (length(find(trigind==ptmp(1)))==0)
                trigind=[trigind;ptmp(1)];

                nt = [nt;labels(ptmp(1))];
                repcnt=[repcnt;repvals(ptmp(end))];
                ptmp = ptmp(1);
                if (strcmp(labels(ptmp),NOTE))
                    toffset = [toffset;tt(ii)*1e3-onsets(ptmp)];
                    ttilend = [ttilend;-tt(ii)*1e3+offsets(ptmp)];
                end
            end
        else
            nfndcnt=nfndcnt+1;
            ptmp = find(offsets>=tt(ii)*1e3);
            [yy,ptmp]=min(abs(0.5*(onsets+offsets)-tt(ii)*1e3));
            if (length(ptmp)>0)
                if (length(find(trigind==ptmp(1)))==0)
                    trigind=[trigind;ptmp(1)];

                    nt = [nt;labels(ptmp(end))];
                    repcnt=[repcnt;repvals(ptmp(end))];
                    ptmp = ptmp(end);
                    if (strcmp(labels(ptmp),NOTE))
                        toffset = [toffset;tt(ii)*1e3-onsets(ptmp)];
                    end
                end
            else
                disp(['BAD Total = ',num2str(nfndcnt)]);
                continue;
            end
        end
        if (length(tdatat)>0)
            ptmp2 = find((tdatat>=onsets(ptmp)*1e-3)&(tdatat<=offsets(ptmp)*1e-3));
            if (length(ptmp2)>0)
                mintmp = [mintmp;min(tdata(ptmp2))];
            end
        end
    end
    
    rd=readrecf(fn);
    trigs(cnt).toffset = toffset;
    trigs(cnt).ttilend = ttilend;
    trigs(cnt).tnotes  = nt;
    trigs(cnt).totalnt = length(labels);
    trigs(cnt).totaltargetnt = length(pos); % total number of target notes
    tmppos = find(onsets*1e-3>=rd.tbefore);
    trigs(cnt).totaltargetnt2 = length(findstr(labels(tmppos),NOTE));
    if (DoLowerCase==1)
        fff = findstr(lower(nt.'),lower(NOTE));
    else
        fff = findstr(nt.',NOTE);
    end
    
    trigs(cnt).matchrcnt=repcnt(fff);
    trigs(cnt).nmatch = length(fff); % number of trigs on target note
    fff2=[1:length(nt)];fff2(fff)=[];
    trigs(cnt).nontrignt=nt(fff2);
    trigs(cnt).nfnd = nfndcnt;
    trigs(cnt).trigmintmp=mintmp; % evtaf trigger valuse min during the note
    nfndcnt=0;
    
    
    if (DoLowerCase==1)
        fff = findstr(lower(labels),lower(NOTE));
    else
        fff = findstr(labels,NOTE);
    end
    ntmintmp=[];
    for ll = 1:length(fff)
        ptmp=find((tdatat>=onsets(fff(ll))*1e-3)&(tdatat<=offsets(fff(ll))*1e-3));
        if (length(ptmp)>0)
            ntmintmp = [ntmintmp;min(tdata(ptmp))];
        end
    end
    trigs(cnt).ntmintmp=ntmintmp;
end
fclose(fid);

vals = zeros([length(trigs),3]);
for ii =1:length(trigs)
    vals(ii,:) = [trigs(ii).nmatch,trigs(ii).totaltargetnt,...
        length(trigs(ii).ttimes)];
end
return;

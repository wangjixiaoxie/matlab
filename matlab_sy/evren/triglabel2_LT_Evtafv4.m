function [vals,trigs]=triglabel(batch,NOTE,PRENT,PSTNT,UseDotTrigFile,DoLowerCase,NoDoCatch,ADDX);


%% LT 7/22/14 - modified to work correctly with birdtaf templates:
% 1) IMPORTANT: KEEP DoLowerCase =0 (i.e. keep low case low and high case high). I did not
% modify code to allow "1";
% 2) If using PSTNT (i.e. equals something other than ''), then look at
% code closely.  should probably work but I don't plan to use that often.


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


searchstr=[PRENT NOTE PSTNT];

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
	    disp(['Hey, .not.mat file has no labels']);
	    continue;
    end

    if (NoDoCatch==1)
	    rdtmp=[];
	    rdtmp=readrecf_LT_evtafv4(fn);
	    if (length(rdtmp)>0)
		    if (rdtmp.iscatch)
			    continue;
		    end
	    end
    end
    
    if (DoLowerCase==1)
        pos = findstr(lower(labels),lower(searchstr));
    else
        pos = findstr(labels,searchstr);
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
        %ppf = findstr(fn,'.ebin');
        %rd=readrecf([fn(1:ppf(end)),'rec_new']);
	if (ADDX==1)
	    ptmp=findstr(fn,'.cbin');
	    if (length(ptmp)==0)
	        ptmp=findstr(fn,'.ebin');
	    end
	    fntemp=[fn(1:ptmp(end)-1),'X.rec'];
	else
	    fntemp=fn;
	end
        rd=readrecf_LT_evtafv4(fntemp);
	if (~isfield(rd,'ttimes'))
		rd.ttimes=[];
	end
        if (strcmp(ext,'.ebin'))
            tt = rd.ttimes/1e3;
        else
            tt = rd.ttimes/1e3;
        end
    end
    fpp = findstr(fn,ext);
    if (ADDX==1)
    	fnt = [fn(1:fpp(end)-1),'X.tmp'];
    else
    	fnt = [fn(1:fpp(end)),'tmp'];
    end
    if (exist(fnt,'file'))
        tdata=load(fnt);
        tdatat = [0:length(tdata)-1]*256/fs + 4.0;
    else
	    tdatat=[];
    end

    
    
    
    %[sm,sp,t,f]=evsmooth(dat,fs,0.01);
    cnt=cnt+1;
   
    repvals = get_repeats(labels);
    
    
    hr=fn2date(fn);
    trigs(cnt).fn=fn;
    trigs(cnt).ftime=hr;
    trigs(cnt).ttimes=tt;
    nt=[];nfndcnt=0;mintmp=[];toffset=[];ttilend=[];repcnt=[];
    trigind=[];
    
    pre_count=length(PRENT); % LT, to use later. takes into account birdtaf.
    post_count=length(PSTNT);
        
    
    
    %loop through all the trigger times
    for ii = 1:length(tt)
        ptmp=find((onsets<=tt(ii)*1e3)&(offsets>=tt(ii)*1e3)); % position of the note (label) that contains the trigger time.
        if (length(ptmp)>0)
            if (length(find(trigind==ptmp(1)))==0)
                trigind=[trigind;ptmp(1)];
                
                nt = [nt;labels(ptmp(1)-pre_count:ptmp(1)+post_count)]; % LT modified
                repcnt=[repcnt;repvals(ptmp(end))];
                ptmp = ptmp(1); % the first instance that it hits?
                % if (strcmp(labels(ptmp),searchstr))
                if (strcmp(labels(ptmp(1)-pre_count:ptmp(1)+post_count),searchstr)); % LT modified 4/12 - previously would not give me offsets if using birdtaf
                    toffset = [toffset;tt(ii)*1e3-onsets(ptmp)]; 
                    ttilend = [ttilend;-tt(ii)*1e3+offsets(ptmp)];
                end
            end
        else
            nfndcnt=nfndcnt+1;
            ptmp = find(offsets>=tt(ii)*1e3);
            [yy,ptmp]=min(abs(0.5*(onsets+offsets)-tt(ii)*1e3)); % takes the average, i.e. the syl that is closest to trigger
            if (length(ptmp)>0)
                if (length(find(trigind==ptmp(1)))==0)
                    trigind=[trigind;ptmp(1)];

                    nt = [nt;labels(ptmp(end)-pre_count:ptmp(end)+post_count)]; % LT modified
                    repcnt=[repcnt;repvals(ptmp(end))];
                    ptmp = ptmp(end);
                    if (strcmp(labels(ptmp(end)-pre_count:ptmp(end)+post_count),searchstr))
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
    
    rd=readrecf_LT_evtafv4(fn);
    trigs(cnt).toffset = toffset;
    trigs(cnt).ttilend = ttilend;
    trigs(cnt).tnotes  = nt;
    trigs(cnt).totalnt = length(labels);
    trigs(cnt).totaltargetnt = length(pos); % total number of target notes
    tmppos = find(onsets*1e-3>=rd.tbefore);
    trigs(cnt).totaltargetnt2 = length(findstr(labels(tmppos),NOTE)); % only the target note, ignoring sequence contingencies
    % LT added 10/31/14 - gives indices of trigs, labels, and 0s and 1s for
    % labels that were hit or escaped.
    trigs(cnt).TrigInds=trigind; 
    trigs(cnt).LabelInds=(pos+length(PRENT))';
    trigs(cnt).Labeled_HitsEscapes = ismember(trigs(cnt).LabelInds, trigs(cnt).TrigInds);
    
    fff=[];
    if (DoLowerCase==1)
        disp('CAUTION: in triglabel2, does not support having (DoLowerCase==1)'); % LT
        fff = findstr(lower(nt.'),lower(NOTE)); % will not work well. see immediately below.
    else
%         fff = findstr(nt.',NOTE); 
        for nn=1:size(nt,1); % LT modified to take into acount sequence
            if strcmp(nt(nn,:),searchstr);
                fff=[fff, nn];
            end
        end
    end
    
    trigs(cnt).matchrcnt=repcnt(fff);
    trigs(cnt).nmatch = length(fff); % number of trigs on target note
    fff2=[1:size(nt,1)];fff2(fff)=[];
    trigs(cnt).nontrignt=nt(fff2,:); % notes that were trigged, but were not desired notes.
    trigs(cnt).nfnd = nfndcnt;
    trigs(cnt).trigmintmp=mintmp; % evtaf trigger valuse min during the note - requires having .tmp file loaded.
    nfndcnt=0;
    
    
    if (DoLowerCase==1)
        disp('CAUTION: in triglabel2, does not support having (DoLowerCase==1)'); % LT
        fff = findstr(lower(labels),lower(NOTE)); % will not work well. see immediately below.
    else
        fff = findstr(labels,searchstr);
        fff = fff + pre_count; % LT: i.e. is there is a preceding note.
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

% close
fclose(fid);

vals = zeros([length(trigs),3]);
for ii =1:length(trigs)
    vals(ii,:) = [trigs(ii).nmatch,trigs(ii).totaltargetnt,...
        length(trigs(ii).ttimes)];
end
return;

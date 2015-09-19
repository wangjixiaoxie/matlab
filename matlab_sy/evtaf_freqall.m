function [fvalsstr]=findwnote4(batch,NOTE,PRENOTE,POSTNOTE,...
                    TIMESHIFT,FVALBND,NFFT,USEFIT,CHANSPEC,ADDX);
% vals=evtaf_freq(bt,fbins,NT,NFFT,CS,USEX);
% vals=[vals;fn2datenum(fn),tempvals,IND,filenum];
%
%  bt - batch file
%  fbins - [Min Freq, Max Freq] to search for peaks (in Hz)
%  NT - target note
%  NFFT - length of template 
%  CS - chan spec
%  USEX - if == 1 look in X.rec for trigger times
% returns vals - 
%   vals=[datenum of file , FREQ vals , Note Index , file number];
% usage example:
%  vals=evtaf_freq('batch.train',[5000,6000],'a',128,'obs0',0);

fvalsstr=[];
if (~exist('PRENOTE'))
    PRENOTE='';
elseif (length(PRENOTE)<1)
    PRENOTE='';
end

if (~exist('CHANSPEC'))
    CHANSPEC='obs0';
elseif (length(CHANSPEC)<1)
    CHANSPEC='obs0';
end

if (~exist('NFFT'))
    NFFT=1024;
elseif (length(NFFT)<1)
    NFFT=1024;
end

if (~exist('USEFIT'))
    USEFIT=1;
elseif (length(USEFIT)<1)
    USEFIT=1;
end

if (~exist('ADDX'))
    ADDX=0;
elseif (length(ADDX)<1)
    ADDX=0;
end

note_cnt=0;avnote_cnt=0;fcnt=0;
ff=load_batchf(batch);
for ifn=1:length(ff)
    fn=ff(ifn).name;
    fnn=[fn,'.not.mat'];
    if (~exist(fnn,'file'))
        continue;
    end
    disp(fn);
    load(fnn)
    labels = lower(labels);
    labels(findstr(labels,'0'))='-';
    if (ADDX==1)
        ptemp=findstr(fnn,'.cbin');
	if (length(ptemp)==0)
        	ptemp=findstr(fnn,'.cbin');
	end
		
        fnrt=[fnn(1:ptemp(end)-1),'X.rec'];
    else
	fnrt=fn;
    end
    rd = readrecf(fnrt);

    [pthstr,tnm,ext] = fileparts(fn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn,CHANSPEC);
    else
        [dat,fs]=evsoundin('',fn,CHANSPEC);
    end

    if (length(dat)==0)
        disp(['hey no data!']);
        continue;
    end
    fcnt=fcnt+1;

    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE])+length(PRENOTE);

    
    %looping through each note in the file
    for ii = 1:length(p)
        ton=onsets(p(ii));toff=offsets(p(ii));

        if (isfield(rd,'ttimes'))
            trigindtmp=find((rd.ttimes>=ton)&(rd.ttimes<=toff));
           
            if (length(trigindtmp)>0)
                TRIG=1;
                if (isfield(rd,'catch'))
                    ISCATCH=rd.catch(trigindtmp);
                else
                    ISCATCH=-1;
                end


            else
                TRIG=0;
                ISCATCH=-1;
            end
        else
            TRIG=0;
            ISCATCH=-1;
        end


%Here frequency is calculated relative to onset of note,but in
%this version I want to calculate frequency relative to detection of
%virtual trigger
        ti1=ceil((TIMESHIFT + ton*1e-3)*fs);
        if (ti1+NFFT-1<=length(dat))
            note_cnt = note_cnt + 1
            dattmp=dat([ti1:(ti1+NFFT-1)]);
            fdattmp=abs(fft(dattmp.*hamming(length(dattmp))));

            %get the freq vals in Hertz
            fvals=[0:length(fdattmp)/2]*fs/(length(fdattmp));
	    fdattmp=fdattmp(1:end/2);
            mxtmpvec=zeros([1,size(FVALBND,1)]);
            
            
            %extra frequency calculation for all these notes
            if(trigindtmp)    
                trigtm=rd.ttimes(trigindtmp)*1e-3;
                inds=fix(trigtm*fs)+[-(NFFT-1):0]'
                dat2=dat(inds);
                fdat=abs(fft(hamming(NFFT).*dat2));
                ffv=get_fft_freqs(NFFT,fs);
                ffv=ffv(1:end/2);
                
                tempvals=[];
                for kk=1:1
                    inds2=find((ffv>=FVALBND(kk,1))&(ffv<=FVALBND(kk,2)));
                    [y,i]=max(fdat(inds2));
                    i=i+inds2(1)-1;
                    i=i+[-3:3]
                    tempvals=[tempvals,sum(ffv(i).*fdat(i).')./sum(fdat(i))];
                end
                fvalsstr(note_cnt).fvirt=tempvals;
            else
                fvalsstr(note_cnt).fvirt=0;
                
                
            end
                %here is where I calculate frequency.

            
            
            
            
            
            for kk = 1:size(FVALBND,1)
                tmpinds=find((fvals>=FVALBND(kk,1))&(fvals<=FVALBND(kk,2)));

		NPNTS=10;
                [tmp,pf] = max(fdattmp(tmpinds));
		pf = pf + tmpinds(1) - 1;
                if (USEFIT==1)
                        tmpxv=pf + [-NPNTS:NPNTS];
			tmpxv=tmpxv(find((tmpxv>0)&(tmpxv<=length(fvals))));

			mxtmpvec(kk)=fvals(tmpxv)*fdattmp(tmpxv);
			mxtmpvec(kk)=mxtmpvec(kk)./sum(fdattmp(tmpxv));
		else
                	mxtmpvec(kk) = fvals(pf);
                end
            end

            [hr,dy,mn,yr]=fn2date(fn);
            fvalsstr(note_cnt).fn     = fn;
            fvalsstr(note_cnt).hour   = hr;
            fvalsstr(note_cnt).day    = dy;
            fvalsstr(note_cnt).month  = mn;
            fvalsstr(note_cnt).mxvals = [1,mxtmpvec];
            fvalsstr(note_cnt).TRIG   = TRIG;
            fvalsstr(note_cnt).CATCH  = ISCATCH;
            fvalsstr(note_cnt).fdat   = fdattmp;
            fvalsstr(note_cnt).datt   = dattmp;
            fvalsstr(note_cnt).ons    = onsets;
            fvalsstr(note_cnt).offs   = offsets;
            fvalsstr(note_cnt).lbl    = labels;
            fvalsstr(note_cnt).ind    = p(ii);
            fvalsstr(note_cnt).ver    = 4;
            fvalsstr(note_cnt).NPNTS  = NPNTS;
        else
            disp('hey');
        end
    end
end
return;
    
    
    
    
    Nbins=3;
NFFT=NFFT*2;
vals=[];
ff=load_batchf(bt);
for ii=1:length(ff)
	fn=ff(ii).name;
	ppp=findstr(fn,'.cbin');
	pppp=findstr(fn,'.');
	tmp=find(pppp<ppp);pppp=pppp(tmp(end));
	filenum=str2num(fn(pppp+1:ppp-1));
	if (~exist(fn,'file'))
		continue;
	end
	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	rd=readrecf(fn,USEX);
	load([fn,'.not.mat']);

	[dat,fs]=evsoundin('',fn,CS);

	if (~isfield(rd,'ttimes'))
		tt=[];
	else
		tt=rd.ttimes*1e-3;
	end
	for jj=1:length(tt)
		pp=find((onsets<=tt(jj)*1e3)&(offsets>=tt(jj)*1e3));
		if (~strcmp(labels(pp),NT))
			continue;
		end
		IND=pp(1);
		inds=fix(tt(jj)*fs)+[-(NFFT-1):0];
		dat2=dat(inds);
		fdat=abs(fft(hamming(NFFT).*dat2));

		ffv=get_fft_freqs(NFFT,fs);
		ffv=ffv(1:end/2);
		
		tempvals=[];
		for kk=1:1%size(fbins,1)
			inds2=find((ffv>=fbins(kk,1))&(ffv<=fbins(kk,2)));
			[y,i]=max(fdat(inds2));
			i=i+inds2(1)-1;
			i=i+[-Nbins:Nbins];
			tempvals=[tempvals,sum(ffv(i).*fdat(i).')./sum(fdat(i))];
		end
		vals=[vals;fn2datenum(fn),tempvals,IND,filenum];
	end
end
return;

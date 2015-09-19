%% LT 7/22/14 - modified:
% 1) to allow choice of evtafv4 or evtaf_amp (only affects getting time info from filename at very end).
% 2) looked through, and all functionality shoudl work.
% 3) uses readrecf_LT_evtafv4, which gets trignotes information as well.
% only imrpotant iuf using EvTAFv4Sim to make rec files.

function vals=evtaf_freq2(bt,fbins,NT,PRENT, PSTNT, NFFT,CS,USEX,evtaf_ver);
%modified 3.4.08 to allow for insertion of PRENT and PSTNT

% vals=evtaf_freq(bt,fbins,NT,NFFT,CS,USEX);
% vals=[vals;fn2datenum(fn),tempvals,IND,filenum];
%
%  bt - batch file
%  fbins - [Min Freq, Max Freq] to search for peaks (in Hz)
%  NT - target note
%  NFFT - length of template 
%  CS - chan spec
%  USEX - if == 1 look in X.rec for trigger times
% evtaf_ver - added by LT, ('amp' for amp; 'v4' for v4)

% returns vals - 
%   vals=[datenum of file , FREQ vals , Note Index , file number];
% usage example:
%  vals=evtaf_freq('batch.train',[5000,6000],'a',128,'obs0',0);


if (~exist('PRENT'))
    PRENT='';
elseif (length(PRENT)<1)
    PRENT='';
end

Nbins=3;
NFFT=NFFT*2;
vals=[];
ff=load_batchf(bt);
for ii=1:length(ff)
% 	ii
    fn=ff(ii).name
	ppp=findstr(fn,'.cbin');
	pppp=findstr(fn,'.');
	tmp=find(pppp<ppp);pppp=pppp(tmp(end));
	filenum=str2num(fn(pppp+1:ppp-1)); % pulls out file number
	if (~exist(fn,'file'))
		continue;
	end
	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	rd=readrecf_LT_evtafv4(fn,USEX);
	load([fn,'.not.mat']);

	[dat,fs]=evsoundin('',fn,CS);

%     labels = lower(labels); LT removed
    labels(findstr(labels,'0'))='-';
	if (~isfield(rd,'ttimes'))
		tt=[];
	else
		tt=rd.ttimes*1e-3;
    end
    p=findstr(labels,[PRENT,NT,PSTNT])+length(PRENT);
    for jj = 1:length(p)  % all labels
        ton=onsets(p(jj));toff=offsets(p(jj));
    
		pp=find((tt*1e3>=ton)&(tt*1e3<=toff));  % picks out the trigger time that fits that labeled syl.
        
        if(~isempty(pp))
            IND=pp(1); % if multiple triggers in one syl, picks out the first one.
            inds=fix(tt(IND)*fs)+[-(NFFT-1):0];  % inds are the 0.008s right before the trigger (trigger is the exact time bin where counters passed threshold)
            dat2=dat(inds);
            fdat=abs(fft(hamming(NFFT).*dat2));

            ffv=get_fft_freqs(NFFT,fs); % frequency bins are 125hz (for NFFT of 128, with 32000 sample rate)
            ffv=ffv(1:end/2);
		
            tempvals=[];
            for kk=1:1%size(fbins,1) % could change to look at different frequency bands.
                inds2=find((ffv>=fbins(kk,1))&(ffv<=fbins(kk,2)));
                [y,i]=max(fdat(inds2)); % picks out max peak
                i=i+inds2(1)-1;
                i=i+[-Nbins:Nbins];
                tempvals=[tempvals,sum(ffv(i).*fdat(i).')./sum(fdat(i))];   % takes weighted sum of frequency values (weighted by power)
            end
            vals=[vals;fn2datenum_eftafv4_lt(fn),tempvals,IND,filenum, p(jj)]; % LT modified to use with either evtaf4 or Amp. (not fn2datenum); p(jj) gives the index fo the label within the song               
        else
            continue;
        end
        
        end
end
return;

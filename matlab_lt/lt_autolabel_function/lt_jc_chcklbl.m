%% LT 10/1/15 - adding description
% INPUTS:
% bt batch file
% NT target note
% PRETIME - duration (sec) before syl
% POSTTIME - duyration post syl to take
% PRENT - pre note, not sure if works perfectly
% POSTNT - "
% CS - leave as '' (is filetype for evsoundin)

% OUTPUTS:
% fnames = filename of .wav file that contains the collected labeled syls
% sylnum = total number of labeled syls found (across all files)

%% LT 10/31/14 - modified to have output wav file include name of syl and timestamp, and to output fnames
% Also gives number of syls found.
function [fnames, sylnum]=lt_jc_chcklbl(bt, NT, PRETIME, POSTTIME, PRENT, POSTNT, CS)
%creates wav file containing all labelled syllables from cbin files
fnames={};

if (~exist('CS'))
    CS='obs0';
elseif (length(CS)==0)
    CS='obs0';
end

if (~exist('PRENT'))
    PRENT='';
elseif (length(PRENT)==0)
    PRENT='';
end

if (~exist('POSTNT'))
    POSTNT='';
elseif (length(POSTNT)==0)
    POSTNT='';
end


%open batch file and get file names
fid=fopen(bt,'r');
files=[];cnt=0;
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    cnt=cnt+1;
    files(cnt).fn=fn;
end
fclose(fid);

%find labelled syllables in each sound file, cut out bit that corresponds
%to label and splice together into one vector
syllwv1=[];
sylnum=0;
for ii = 1:length(files);
    fn=files(ii).fn;
    if exist([fn,'.not.mat'],'file')
        load([fn,'.not.mat']);
    else
        labels=[];
    end
    labels(findstr(labels,'0'))='-';
    
    
    %         if (ii==1)
    %         [dat,fs]=evsoundin('',fn,CS);
    %
    %     end
    
    [dat,fs]=evsoundin('',fn,CS);
    dat = dat *1e-4;
    
    NPRE=ceil(PRETIME*fs);%time before syllable in samples
    NPOST=ceil(POSTTIME*fs);%time after syllable in samples
    srchstr=[PRENT,NT,POSTNT];

    pp=findstr(labels,srchstr)+length(PRENT);
    if ~isempty(pp)>0;
        for jj=1:length(pp)
            onind=fix(round(onsets(pp(jj))*1e-3*fs));%onset time into samples
            enind = fix(round(offsets(pp(jj))*1e-3*fs));%offset time into samples
            st=onind;%onset time - pretime in samples
            en=enind;%offset time + posttime in samples
            
            if (st<1)
                st=1;
            end
            if (en>length(dat))
                en=length(dat);
            end
            
            dat_tmp=dat(st:en);
            dat_tmp = [zeros(NPRE,1); dat_tmp; zeros(NPOST,1)];
            syllwv1 = cat(1,syllwv1, dat_tmp);
            sylnum=sylnum+1;
        end
    end
end


tstamp=lt_get_timestamp(0);
fnames{1}=['syllwv_' NT '_' tstamp '.wav'];
% fnames{2}=['syllwv_' NT '_' tstamp 'B.wav'];
wavwrite(syllwv1,fs,fnames{1});
% wavwrite(syllwv2, fs,fnames{2});



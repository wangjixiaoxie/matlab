%% LT 7/22/14 - modified to accomodate evtafv4 column syntax:
% e.g. (a+b+c)*(d+e+f) in Evtafv4 means (a or b or c) and (d or e or f).  That contrasts with evtaf_amp,
% where a or b or c and d or e or f meant {(a or b or c) and (d)} or e or f (i.e. any time an or occurs it lumps everything to the left of it as a single entity)
% (causing a limitation where at least one side of an "and" must be a single column.

% triglabel originally goes left to right.  I will modify so that each
% bracketed entity is performed first, then their outputs are compared.

% NOTE: only supports one layer of recursiveness - e.g. (a+b)*(c+d) and not
% ((a+b)*(c+d))+f
% ( i'm not even sure if evtafv4 itself supports that

% NOTE: added arg col_logic below

% NOTE: Should have all stuff within brackets have same logic (e.g. not
% (a+b*c).  I am not sure what the order of operations is in actual evtaf,
% so not sure how to code here. regardless that is not necessary as can
% just put brackets around things to be unambiguous (e.g. (a+b)*c or
% a+(b*c));


function trigs=get_trigt(bt,cntrng,refrac,NFFT,ADDX,TEMPFILE,col_logic);
%trigs=get_trigt2(bt,cntrng,refrac,NFFT,ADDX,TEMPFILE);
%
% bt=batch file
% cntrng = struct with MIN, MAX, TH, NOT, BTMIN, AND (==1 means do and)
% refrac in seconds
% NFFT is the number of points in the templ file
%ADDX if ==1 add the X to the tmp file name to use
% termperary tmp files defualt == 1
% TEMPFILE make a .###X.rec file to avoid overwrite
%    if this is == 1 (defualt)
% DO_OR if ==1 then will use ORs between templates
% col_logic='(a+b+c)*(d+f)' --> as a string, no spaces.  identical to what used in evtafv4


% ADDED THE NEW COUNTER VERSION WERE IS THE COUNTER HAS
% TO BE >=N and val<TH for the count to continue goin up


if (~exist('ADDX'))
    ADDX=1;
end
if (~exist('TEMPFILE'))
    TEMPFILE=1;
end
if (~exist('NFFT'))
    NFFT=128;
end

if (~exist('refrac'))
    refrac=100e-3;
end

%if (~exist('DO_OR'))
%    DO_OR=0;
%else
%    if (length(DO_OR)==0)
%        DO_OR=0;
%    end
%end

for ii=1:length(cntrng)
    if (~isfield(cntrng(ii),'NOT'))
        cntrng(ii).NOT=0;
    end
    if (~isfield(cntrng(ii),'MODE'))
        cntrng(ii).MODE=1;
    end
end

files=load_batchf(bt);

NCNT=length(cntrng);
cnt=zeros([1,NCNT]);
trigs=[];

%% extracting logic in usable form from col_logic

% find the locations of the bracketed groups
open_brackets=findstr(col_logic,'(');
close_brackets=findstr(col_logic,')');

num_groups=size(open_brackets,2); % number of groups within different paranthese

for iii=1:num_groups;
    group_contents{iii}=col_logic(open_brackets(iii)+1:close_brackets(iii)-1); % what are string contents of each group
    group_num_columns{iii}=(size(group_contents{iii},2)+1)/2; % how many columns per group?
    % get and/or logics between columns within groups
    for jjj=1:group_num_columns{iii}-1; 
        if strcmp(group_contents{iii}(jjj*2),'+'); % + means or 
            and_within_group{iii}(jjj)=0; % i.e. 'or'
        elseif strcmp(group_contents{iii}(jjj*2),'*'); 
            and_within_group{iii}(jjj)=1; % i.e. 'and'
        end
    end
    
    % within groups, convert from letters to numbers corresponding to
    % columns position.
    for jjj=1:group_num_columns{iii};
        columns{iii}(jjj)=double(group_contents{iii}(jjj*2-1))-double('a')+1; % converts a to 1, b to 2, etc
    end
    
    % get and/or logic between groups.
    if iii<num_groups;
        if strcmp(col_logic(close_brackets(iii)+1),'+');
            and_between_groups(iii)=0;
        elseif strcmp(col_logic(close_brackets(iii)+1),'*');
            and_between_groups(iii)=1;
        end
    end
end

% OUTPUT:
% columns; 
% and_within_group;
% and_between_groups;

%% Get trigs

for ijk=1:length(files)
    fn=files(ijk).name;
    if (~exist(fn,'file'))
        continue;
    end
    
    pp=findstr(fn,'.cbin');
    if (length(pp)==0)
        pp=findstr(fn,'.ebin');
    end
    
    if (ADDX==1)
        fnt=[fn(1:pp(end)-1),'X.tmp'];
    else
        fnt=[fn(1:pp(end)),'tmp'];
    end
    
    if (~exist(fnt,'file'))
        disp('hey, X.tmp does not exist - should run mk_tempf to make those');
        continue;
    end
    disp(fnt);
    
    rdat=readrecf(fn);
    tmpdat=load(fnt); % loads temp file
    tmpdat2=zeros([fix(length(tmpdat)/NCNT),NCNT]);
    for ii=1:NCNT % converts temp file to multi-column format
        placeholder=tmpdat(ii:NCNT:end);
        tmpdat2(:,ii)=placeholder(1:size(tmpdat2,1));
    end
    tmpdat=tmpdat2;
    clear tmpdat2;
    
    fs=rdat.adfreq;
    %TH=1.5;fs=32000;
    refracsam=ceil(refrac/(2*NFFT/fs));  % how many 0.008s samples are there in refract
    lasttrig=-refrac;tt=[];
    cnt=0*cnt;
    for kk=1:NCNT
        if (cntrng(kk).MODE==0)
            cnt(kk)=cntrng(kk).MAX+1;
        end
    end
    
    cntvals=zeros([NCNT,size(tmpdat,1)]);
    for ii = 1:size(tmpdat,1)
        for kk=1:NCNT
            if (tmpdat(ii,kk)<=cntrng(kk).TH)
                if (cntrng(kk).MODE==1)
                    cnt(kk)=cnt(kk)+1;
                else
                    if (cnt(kk)>=cntrng(kk).BTMIN)
                        cnt(kk)=0;
                    else
                        cnt(kk)=cnt(kk)+1;
                    end
                end
            else
                if (cntrng(kk).MODE==0)
                    cnt(kk)=cnt(kk)+1;
                else
                    cnt(kk)=0;
                end
            end
        end
        cntvals(:,ii)=cnt.';
        
        %if (DO_OR==0)
        %    trig=1;
        %else
        %    trig=0;
        %end
        
        %% LT get trigs using evtafv4 logic (MISSING REPEATS)
        for iii=1:num_groups;
            for kkk=1:group_num_columns{iii}; % number of cols in group
                kk=columns{iii}(kkk); % kk= original column index.
                if ((cnt(kk)>=cntrng(kk).MIN)&(cnt(kk)<=cntrng(kk).MAX));
                    ntrig=1; % satisfy threshold
                else
                    ntrig=0;
                end
                if (cntrng(kk).NOT==1)
                    ntrig=~ntrig;
                end
                
                if (kkk==1) % if first, then no previous logic to use
                    trig=ntrig;
                else
                    if and_within_group{iii}(kkk-1)==1; 
                        trig = trig & ntrig; 
                    else
                        trig = trig | ntrig;
                    end
                end
            end
            trig_per_group(iii)=trig; % output for each group (1 means trig, 0 no)
        end
        
        % combine groups using logic between groups
        for iii=1:num_groups;
            ntrig=trig_per_group(iii); 
            
            if iii==1;
                trig=ntrig;
            else
                if and_between_groups(iii-1)==1;
                    trig=trig & ntrig;
                else
                    trig=trig | ntrig; % final output, trig is either 0 or 1 (hit), for a given time bin in the tmp file.
                end
            end
        end
        
        
        %%
        if (trig) % output time is the time at the right edge of the trigger bin.
            if (abs(ii-lasttrig)>refracsam)
                tt=[tt;((ii*NFFT*2/fs)+rdat.tbefore)*1e3];
                lasttrig=ii;
            end
        end
    end
    trigs(ijk).fn=fn;
    trigs(ijk).cntvals=cntvals;
    
    rdat.ttimes=tt;
    tmp=zeros([length(cntrng),1]);
    for ii=1:length(tmp)
        tmp(ii)=cntrng(ii).TH;
    end
    rdat.thresh=tmp;
    if (TEMPFILE==1)
        fntemp=fn;
        pp=findstr(fn,'.cbin');
        if (length(pp)==0)
            pp=findstr(fn,'.bbin');
        end
        if (length(pp)==0)
            pp=findstr(fn,'.ebin');
        end
        if (length(pp)>0)
            fntemp=[fn(1:pp(end)-1),'X',fn(pp(end):end)];
            wrtrecf(fntemp,rdat);
        end
    else
        wrtrecf(fn,rdat);
    end
end
return
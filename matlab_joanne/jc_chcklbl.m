function jc_chcklbl(bt, NT, PRETIME, POSTTIME, PRENT, POSTNT, CS)
%creates wav file containing all labelled syllables from cbin files

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
syllwv1 = [];
syllwv2 = [];
if mod(length(files),2) == 0
    
    for ii=1:length(files)/2
        fn=files(ii).fn;
        if strfind(fn,'filt') > 0
            load([fn(1:end-5),'.not.mat'])

            elseif (exist([fn,'.not.mat'],'file'))
                load([fn,'.not.mat']);
        else
            labels=[];
        end
        labels(findstr(labels,'0'))='-';

        if (ii==1)
            [dat,fs]=evsoundin('',fn,CS);

            NPRE=ceil(PRETIME*fs);%time before syllable in samples
            NPOST=ceil(POSTTIME*fs);%time after syllable in samples
        end

%         disp(fn);
        srchstr=[PRENT,NT,POSTNT];
        pp=findstr(labels,srchstr)+length(PRENT);
        if (length(pp)>0)

            [filepath,filename,fileext] = fileparts(fn);
            if(strcmpi(fileext,'.wav'))
                [dat,fs] = wavread(fn);
                dat = dat *10e3; % boost amplitude to cbin levels
            else
                [dat,fs]=evsoundin('',fn,CS);
                dat = dat *1e-4;
            end
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
            end
        end
    end
    
    for ii = (length(files)/2)+1:length(files)
        fn=files(ii).fn;
        if strfind(fn,'filt') > 0
            load([fn(1:end-5),'.not.mat'])

            elseif (exist([fn,'.not.mat'],'file'))
                load([fn,'.not.mat']);
        else
            labels=[];
        end
        labels(findstr(labels,'0'))='-';

        if (ii==1)
            [dat,fs]=evsoundin('',fn,CS);

            NPRE=ceil(PRETIME*fs);%time before syllable in samples
            NPOST=ceil(POSTTIME*fs);%time after syllable in samples
        end

%         disp(fn);
        srchstr=[PRENT,NT,POSTNT];
        pp=findstr(labels,srchstr)+length(PRENT);
        if (length(pp)>0)

            [filepath,filename,fileext] = fileparts(fn);
            if(strcmpi(fileext,'.wav'))
                [dat,fs] = wavread(fn);
                dat = dat *10e3; % boost amplitude to cbin levels
            else
                [dat,fs]=evsoundin('',fn,CS);
                dat = dat *1e-4;
            end
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

                syllwv2 = cat(1,syllwv2, dat_tmp);
            end
        end
    end
    
else
    for ii = 1:floor(length(files)/2)
        fn=files(ii).fn;
        if strfind(fn,'filt') > 0
            load([fn(1:end-5),'.not.mat'])

            elseif (exist([fn,'.not.mat'],'file'))
                load([fn,'.not.mat']);
        else
            labels=[];
        end
        labels(findstr(labels,'0'))='-';

        if (ii==1)
            [dat,fs]=evsoundin('',fn,CS);

            NPRE=ceil(PRETIME*fs);%time before syllable in samples
            NPOST=ceil(POSTTIME*fs);%time after syllable in samples
        end

%         disp(fn);
        srchstr=[PRENT,NT,POSTNT];
        pp=findstr(labels,srchstr)+length(PRENT);
        if (length(pp)>0)

            [filepath,filename,fileext] = fileparts(fn);
            if(strcmpi(fileext,'.wav'))
                [dat,fs] = wavread(fn);
                dat = dat *10e3; % boost amplitude to cbin levels
            else
                [dat,fs]=evsoundin('',fn,CS);
                dat = dat *1e-4;
            end
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
            end
        end
    end
    
    for ii = ceil(length(files)/2):length(files)
        fn=files(ii).fn;
        if strfind(fn,'filt') > 0
            load([fn(1:end-5),'.not.mat'])

            elseif (exist([fn,'.not.mat'],'file'))
                load([fn,'.not.mat']);
        else
            labels=[];
        end
        labels(findstr(labels,'0'))='-';

        if (ii==1)
            [dat,fs]=evsoundin('',fn,CS);

            NPRE=ceil(PRETIME*fs);%time before syllable in samples
            NPOST=ceil(POSTTIME*fs);%time after syllable in samples
        end

%         disp(fn);
        srchstr=[PRENT,NT,POSTNT];
        pp=findstr(labels,srchstr)+length(PRENT);
        if (length(pp)>0)

            [filepath,filename,fileext] = fileparts(fn);
            if(strcmpi(fileext,'.wav'))
                [dat,fs] = wavread(fn);
                dat = dat *10e3; % boost amplitude to cbin levels
            else
                [dat,fs]=evsoundin('',fn,CS);
                dat = dat *1e-4;
            end
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

                syllwv2 = cat(1,syllwv2, dat_tmp);
            end
        end
    end
end


wavwrite(syllwv1,fs,'syllwv1.wav');
wavwrite(syllwv2, fs, 'syllwv2.wav');


            
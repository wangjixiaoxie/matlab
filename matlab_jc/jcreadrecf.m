function recdata=readrecf(fname,ADDX);
% recdata=readrecf(fname,ADDX);
%
%     adfreq: 32000
%      nchan: 2
%      nsamp: 199684
%    iscatch: 0
%     ttimes: []
%    tbefore:
%     tafter:
%     thresh:
%    outfile:
%     header:

if (~exist('ADDX'))
	ADDX=0;
else
	if (length(ADDX)==0)
		ADDX=0;
	end
end

pp = findstr(fname,'.rec');
if (length(pp)<1)
    pp3=findstr(fname,'.not.mat');
    if (length(pp3)>0)
        fname=fname(1:pp3(end)-1);    
    end
    pp2 = findstr(fname,'.');
    if (length(pp2)<1)
        recf = [fname,'.rec'];
    else
        recf = [fname(1:pp2(end)),'rec'];
    end
else
    recf = fname;
end

if (ADDX==1)
	pp=findstr(recf,'.rec');
	recf=[recf(1:pp(1)-1),'X.rec'];
end

if (~exist(recf,'file'))
    %disp(['Rec file : ',recf,' does not exist']);
    recdata = [];
    return;
end

fid = fopen(recf,'r');
flsv=[];
while (1)
    if (length(flsv)==0)
        fl = fgetl(fid);
    else
        fl = flsv;
        flsv=[];
    end
    if (~ischar(fl))
        break;
    end
    
    if length(fl)==0
        continue;
    end

    if (length(findstr(upper(fl),'CATCH'))>0)
        pp = findstr(fl,'=');
        recdata.iscatch = str2num(fl(pp(end)+1:end));

    elseif (length(findstr(upper(fl),'CHANS'))>0)
        pp = findstr(fl,'=');
        recdata.nchan = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'ADFREQ'))>0)
        pp = findstr(fl,'=');
        recdata.adfreq = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'SAMPLES'))>0)
        pp = findstr(fl,'=');
        recdata.nsamp = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'T AFTER'))>0)
        pp = findstr(fl,'=');
        recdata.tafter = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'T BEFORE'))>0)
        pp = findstr(fl,'=');
        recdata.tbefore = str2num(fl(pp(end)+1:end));
    elseif (length(findstr(upper(fl),'OUTPUT SOUND FILE'))>0)
        pp = findstr(fl,'=');
        recdata.outfile = fl(pp(end)+1:end);
     elseif (length(findstr(upper(fl),'STIMULUS'))>0)
        pp = findstr(fl,':');
        recdata.outfile = fl(pp(end)+2:end);
    elseif (length(findstr(upper(fl),'THRESHOLDS'))>0)
        tmpvec = [];
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(str2num(fl))==0)
                flsv = fl;
                break;
            end
            tmpvec = [tmpvec;str2num(fl)];
        end
        recdata.thresh = tmpvec;
        if (~ischar(fl))
            break;
        end
    elseif (length(findstr(lower(fl),'feedback information'))>0)
        tmpvec = [];tmpvec2 = [];
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(fl)<2)
                continue;
            end
            tpos = findstr(fl,'msec');
            tmpvec = [tmpvec;str2num(fl(1:tpos(1)-1))];
	    tpos = findstr(fl,':');
	    tmpvec2=[tmpvec2;(length(findstr(lower(fl(tpos+1:end)),'catch'))>0)];
        end
        recdata.ttimes = tmpvec;
        recdata.catch  = tmpvec2;
        if (~ischar(fl))
            break;
        end
    elseif (length(findstr(lower(fl),'trigger times'))>0)
        tmpvec = [];
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(fl)<2)
                continue;
            end
            tmpvec = [tmpvec;str2num(fl)];
        end
        recdata.ttimes = tmpvec*1e3;
        if (~ischar(fl))
            break;
        end
    elseif (length(findstr(upper(fl),'FILE CREATED:'))>0)
        % save the next 4 lines and this line as the header
        headr{1} = fl;

        for ijk=1:4
            fl=fgetl(fid);
            headr{ijk+1}=fl;
        end
        recdata.header=headr;
    end
end

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
timesource=char(headr(1));
pp0=findstr(timesource,'2008');
rawtimestr=timesource(ppp+6:end);
% calculate the time of the song's end (in seconds after midnight)
rawtimeend=str2num(rawtimestr(1:2))*3600+str2num(rawtimestr(4:5))*60+str2num(rawtimestr(7:8)); 
% determine the length of the song in seconds
longsong=char(headr(5));
pp1=findstr(longsong,'=');
pp2=findstr(longsong,'ms');
rawlongstr=longsong(pp1+2:pp2-2);
% determine the time of the song's starting
rawtimestart=rawtimeend-str2num(rawlongstr)/1000;





if exist('recdata')
    if (~isfield(recdata,'ttimes'))
	recdata.ttimes=[];
    end
else
    recdata=[];
    return;
end
    fclose(fid);
return;

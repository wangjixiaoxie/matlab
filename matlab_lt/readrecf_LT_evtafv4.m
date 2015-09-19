%% LT 12/7/14 - modified to read trignotes
% This is now compatible with triglabel2_LT_Evtafv4
% extracts not just trigger times, but also the associated template notes
% (i.e. 0, 1, 2, ...). evtafv4 can use multiple templates.


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

if (~exist('ADDX')) % if on, adds an 'X' before .rec to the file.
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
    elseif (length(findstr(upper(fl),'STIMULUS:'))>0)
        pp = findstr(fl,':');
        recdata.outfile = fl(pp(end)+2:end);
    elseif (length(findstr(upper(fl),'OUTPUT SOUND FILE'))>0)
        pp = findstr(fl,'=');
        recdata.outfile = fl(pp(end)+1:end);
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
        tmpvec = [];tmpvec2 = [];tmpvec3=[];pbname={};
        while (1)
            fl = fgetl(fid);
            if (~ischar(fl))
                break;
            end
            if (length(fl)<2) % if is filler(''), then go to next line, but it is still ttimes.
                continue;
            end
            tpos = findstr(fl,'msec');
            tmpvec = [tmpvec;str2num(fl(1:tpos(1)-1))];
            tpos  = findstr(fl,':');
            tpos2 = findstr(lower(fl(tpos+1:end)),'catch');
            tmpvec2=[tmpvec2;(length(tpos2)>0)];
            
            % get information on note num
            xx=findstr(lower(fl),'templ');
            tmpvec3=[tmpvec3; str2num(fl(xx+8))];
            
            
            if (length(tpos2)>0)
                pbname{length(pbname)+1}=fl(tpos+2:end);
            else
                pbname{length(pbname)+1}=fl(tpos+2:end);
            end
        end
        
        recdata.ttimes = tmpvec;
        recdata.catch  = tmpvec2;
        recdata.trignote = tmpvec3;
        recdata.pbname = pbname;
        
        
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
        
        % LT - used to get trignote, but easier method above, as combined
        % with feeback info.
        
        %     elseif (length(findstr(lower(fl),'template notenum'))>0)
        %         tmpvec = [];
        %         while (1)
        %             fl = fgetl(fid);
        %             if (~ischar(fl))
        %                 break;
        %             end
        %             if (length(fl)<2)
        %                 continue;
        %             end
        %             tmpvec = [tmpvec;str2num(fl)];
        %         end
        %         recdata.trignote = tmpvec;
        %         if (~ischar(fl))
        %             break;
        %         end
        %     end
    end
    if (~isfield(recdata,'ttimes'))
        recdata.ttimes=[];
    end
end

if (~isfield(recdata,'trignote'))
    recdata.trignote=[];
end

fclose(fid);
return;

function [datsvnt,datsv]=repanal(bt,NT);
%[datsvnt,datsv]=repanal(bt,NT);
% datsvnt - data for single specified note
%     [which rep indx in song, fix(label), datenum of song,rep cnt]
%
% datsv - data for all repeated notes
%       [datenum,btfulenum,fix(label), repcnt];
%

trigs=[];
cnt=0;
fid=fopen(bt,'r');
rpnums=[];
datsv=[];
datsvnt=[];filecnt=0;
while (1)
    fn=fgetl(fid);
    if (~ischar(fn));
        break;
    end
    if (~exist([fn,'.not.mat']))
        continue;
    end
    load([fn,'.not.mat']);
    if (length(labels)==0)
        disp(['hey no labels :',fn]);
        continue;
    end

    [h,d,m,y]=fn2date(fn);
    mnt = (h-fix(h))*60;
    dtnm=datenum([y,m,d,h,mnt,0]);

    df = diff(fix(labels));
    pp = find(df~=0);
    filecnt=filecnt+1;
    if (length(pp)==0)
        vals=[1,length(labels),fix(labels(1))];
    else
        vals=zeros([length(pp)+1,3]);
        vals(1,:)=[1,pp(1),fix(labels(1))];
        if (length(pp)>1)
            for ii = 2:length(vals)-1
                vals(ii,:)=[vals(ii-1,2)+1,pp(ii),fix(labels(pp(ii)))];
            end
        end
        vals(end,:)=[vals(end-1,2)+1,length(labels),fix(labels(end))];
    end
    tmpdt=zeros([size(vals(:,3))]);
    datsv=[datsv;tmpdt+dtnm,tmpdt+filecnt,vals(:,3),vals(:,2)-vals(:,1)+1];
    pppp=find(vals(:,3)==fix(NT(1)));
    if (length(pppp)>0)
        datsvnt=[datsvnt; ...
            [1:length(pppp)]',vals(pppp,3),dtnm+zeros([length(pppp),1]),vals(pppp,2)-vals(pppp,1)+1];
    end


    rd=readrecf(fn);
    tt=rd.ttimes;

    tmpvals=[];
    note=[];
    cnt=cnt+1;
    for  ii = 1:length(tt)
        ppp=find((onsets<=tt(ii))&(offsets>=tt(ii)));
        if (length(ppp)==1)
            ppp2=find((vals(:,1)<=ppp)&(vals(:,2)>=ppp));
            if (length(ppp2)==1)
                rep_num = ppp-vals(ppp2,1)+1;
                tmpvals=[tmpvals;rep_num];
                rpnums=[rpnums;rep_num];
                note=[note,vals(ppp2,3)];
            else
                disp(['hey2 ',num2str(cnt)]);
            end
        else
            disp(['hey ',num2str(cnt),' ',fn]);
        end
    end
    trigs(cnt).rpnums=tmpvals;
    trigs(cnt).fn=fn;
    trigs(cnt).tt=tt;
    trigs(cnt).note=char(note);
end
fclose(fid);
return;

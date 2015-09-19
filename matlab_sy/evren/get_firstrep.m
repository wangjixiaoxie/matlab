function [valsv]=get_firstrep(bt);
%[valsv]=get_firstrep(bt);
% valsv -- [fix(label), num repsi, datenum];
%
valsv=[];
fid=fopen(bt,'r');
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    load([fn,'.not.mat']);

    lblnum=fix(labels);

    df = diff(fix(lblnum));
    pp = find(df~=0);
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

    reps = [vals(:,3),vals(:,2)-vals(:,1)+1];
    pp=find((reps(:,1)~=fix('-'))&(reps(:,1)~=fix('0')));
    [tmp,tmp2,tmp3,tmp4,dn]=fn2date(fn);
    if (length(pp)>0)
        valsv=[valsv;reps(pp(1),:),dn];
    end
end
fclose(fid);
return;
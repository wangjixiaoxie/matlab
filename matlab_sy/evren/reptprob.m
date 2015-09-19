function [prob,nt,cnts]=reptprob(bt);
%[prob,nt]=reptprob(bt);
%

cnt=0;
fid=fopen(bt,'r');
nt=[];
while (1)
    fn=fgetl(fid);
    if (~ischar(fn));
        break;
    end
    load([fn,'.not.mat']);
    if (length(labels)==0)
        disp(['hey no labels :',fn]);
        continue;
    end
    cnt=cnt+1;

    [h,d,m,y]=fn2date(fn);
    mnt = (h-fix(h))*60;
    dtnm=datenum([y,m,d,h,mnt,0]);

    df = diff(fix(labels));
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

    % get rid of non-labeled notes
    pp=find((vals(:,3)~=fix('-'))&(vals(:,3)~=fix('0')));
    vals=vals(pp,:);

    tmpvals{cnt}=vals(:,3).';
    
    for ii=1:length(vals(:,3))
        if (length(find(nt==vals(ii,3)))==0)
            nt=[nt,vals(ii,3)];
        end
    end
end
fclose(fid);

nt=sort(nt);
prob=zeros(length(nt));
for ii = 1:length(tmpvals)
    tmpvec=tmpvals{ii};
    for jj = 1:length(tmpvec)-1
        ind1=find(nt==tmpvec(jj));
        ind2=find(nt==tmpvec(jj+1));
        prob(ind1,ind2)=prob(ind1,ind2)+1;
    end
end
cnts=prob;
for ii = 1:length(prob)
    prob(ii,:)=prob(ii,:)./sum(prob(ii,:));
end

return;
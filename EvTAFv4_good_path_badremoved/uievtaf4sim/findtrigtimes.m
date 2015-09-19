function trigtimes=findtrigtimes(hObject,handles);
%

vals=handles.TAFVALS;
ntempl=size(vals,2);

CntStr=handles.CntRng;
tt=[];

%initialize counters
cntr=zeros([1,ntempl]);
for jj=1:ntempl
    if (CntStr(jj).Mode==1)
        cntr(jj)=0;
    else
        cntr(jj)=CntStr(jj).MaxCnt+1;
    end
end

lasttrig = -handles.REFRAC;
for ii = 1:size(vals,1)
    for jj = 1:ntempl
        cval=vals(ii,jj);
        if (cval<=CntStr(jj).Thresh)
            if (CntStr(jj).Mode==0)
                cntr(jj)=0;
            else
                cntr(jj)=cntr(jj)+1;
            end
        else
            if (CntStr(jj).Mode==0)
                cntr(jj)=cntr(jj)+1;
            else
                cntr(jj)=0;
            end
        end
    end

    trignowsv=zeros([1,ntempl]);
    for jj = 1:ntempl
        tmp=((cntr(jj)>=CntStr(jj).MinCnt)&(cntr(jj)<CntStr(jj).MaxCnt));
        if (CntStr(jj).DoNot==1)
            tmp=~tmp;
        end
        trignowsv(jj) = tmp;
    end

    if (ntempl>1)
        if (CntStr(1).DoAND==1)
            trignow = trignowsv(1) & trignowsv(2);
        else
            trignow = trignowsv(1) | trignowsv(2);
        end
        for jj=2:ntempl-1
            if (CntStr(jj).DoAND==1)
                trignow = trignow & trignowsv(jj+1);
            else
                trignow = trignow | trignowsv(jj+1);
            end
        end
    else
        trignow=trignowsv;
    end

    tnow = ii*(2*size(handles.Templ,1)/handles.fs) + handles.TBEFORE;
    if (trignow==1)&((tnow-lasttrig)>=handles.REFRAC)
        lasttrig=tnow;
        tt=[tt;tnow];
    end
end

trigtimes = tt;%*2*size(handles.Templ,1)/handles.fs+handles.TBEFORE;
return
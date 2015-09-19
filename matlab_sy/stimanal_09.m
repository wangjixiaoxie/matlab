clear indcatch indnotcatch
for ii=1:length(fv)
 fvcur=fv{ii}   
 indcatch{ii}=[];
        indnotcatch{ii}=[];   
 for jj=1:length(fvcur)
        
        if(fvcur(jj).CATCH)
            indcatch{ii}=[indcatch{ii} jj]
        else
            indnotcatch{ii}=[indnotcatch{ii} jj]
        end
    end
end

%plot cvs.
figure
for ii=1:4
   ax(ii)=subplot(4,1,ii);
end

for ii=1:8
    stdv_mat{ii}=std(ct{ii},0,2);
    mn_mat{ii}=mean(ct{ii},2);
    cv_mat{ii}=stdv_mat{ii}./mn_mat{ii};
end

for ii=1:8
    axvl=ceil(ii/2);
    if(mod(ii,2))
        col='k'
    else
        col='r'
    end
    axes(ax(axvl))
    plot(pitchtms,mn_mat{ii},'Color',col);
    hold on;
end


%get amplitude matrices.
for ii=1:6
    
    pathvl=avls.pvls{ii}
    cmd=['cd ' pathvl]
    eval(cmd);
    bt=avls.cvl{ii};
[avna,t,f,catchsm,sm_mat{ii}]=get_avn3(bt,'a',0.2,0.2,'','','obs0');
end
figure
for ii=1:4
   ax2(ii)=subplot(4,1,ii);
end
for ii=1:8
    stdv_sm{ii}=std(sm{ii},0,2);
    mn_sm{ii}=mean(sm{ii},2);
    cv_sm{ii}=stdv_sm{ii}./mn_sm{ii};
end
for ii=1:8
    axvl=ceil(ii/2);
    if(mod(ii,2))
        col='k'
    else
        col='r'
    end
    axes(ax2(axvl))
    tms=0:.4/12801:.4-.4/25601
    plot(tms,cv_sm{ii},'Color',col);
    hold on;
end

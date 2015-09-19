function PlotTafVals(hObject,handles);
%

axes(handles.ValAxes);hold off;

fn=handles.files{handles.Nfile};
fs=handles.fs;

[pth,nm,ext]=fileparts(fn);
if (exist([pth,nm,'.tmp'],'file')&(get(handles.UseSim,'Value')==get(handles.UseSim,'Min')))
    %set(handles.FromTempFileBox,'Value',get(handles.FromTempFileBox,'Max'));
    rdata=readrecf(fn);
    
    if (handles.USEXTMP==0)
        tmpdat=load([pth,nm,'.tmp']);
    else
        tmpdat=load([pth,nm,'X.tmp']);
    end
    vals=zeros([fix(length(tmpdat)/handles.NTempl),handles.NTempl]);
    for ii = 1:handles.NTempl;
        vals(:,ii)=tmpdat(ii:(handles.NTempl):end);
    end

    v=axis;
    axes(handles.ValAxes);hold off;
    tmpplt=plot([1:length(vals)].*(2*handles.TemplLen)./fs  + rdata.tbefore,vals,'s-');hold on;
    ylim(v(3:4));
    hold on;
    handles.TBEFORE=rdata.tbefore;
else
    %set(handles.FromTempFileBox,'Value',get(handles.FromTempFileBox,'Min'));
    
    vals=evtafsim(handles.dat,fs,handles.Templ);
    
    v=axis;
    axes(handles.ValAxes);hold off;
    tmpplt=plot([1:length(vals)].*(2*handles.TemplLen)./fs,vals,'s-');
    ylim(v(3:4));hold on;
    handles.TBEFORE=0.0;
end

handles.ValsPltHndl = tmpplt;
t=handles.t;
tmpplt2=zeros([length(handles.CntRng),1]);
for ii = 1:length(handles.CntRng)
    clr=get(tmpplt(ii),'Color');
    tmpplt2(ii)=plot([t(1),t(end)],[1,1]*handles.CntRng(ii).Thresh);hold on;
    set(tmpplt2(ii),'Color',get(tmpplt(ii),'Color'));
end
handles.ThreshLineHand = tmpplt2;
handles.TAFVALS=vals;
guidata(hObject,handles);

axes(handles.LabelAxes);
trigtimes = findtrigtimes(hObject,handles);
handles.TrigT=trigtimes;
guidata(hObject,handles);

%if (handles.TRIGPLT~=-1)
%	delete(handles.TRIGPLT);
%end
%tmp_plt=plot(trigtimes,0*trigtimes-0.5,'bs');
%axis([0,t(end),-1,1]);
%handles.TRIGPLT=tmp_plt;
%guidata(hObject,handles);
return;

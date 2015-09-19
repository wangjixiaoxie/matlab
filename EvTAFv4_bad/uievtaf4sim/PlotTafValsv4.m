function PlotTafValsv4(hObject,handles,REPLOT);
%

axes(handles.ValAxes);hold off;

fn=handles.files{handles.Nfile};
fs=handles.fs;

[pth,nm,ext]=fileparts(fn);

%set(handles.FromTempFileBox,'Value',get(handles.FromTempFileBox,'Min'));
rdata=readrecf(fn);
handles.TBEFORE=handles.OP.FileBufferLeng;
%%%%%%%%%%%%%%%%%%%

if (REPLOT==1)
    [ND,TemplMatchVals]=EvTAFv4(handles.ND,handles.OP,handles.dat,fs);
    handles.ND=ND;
    handles.TAFVALSALL=TemplMatchVals;
    guidata(hObject,handles);
    temptrigs=[];
    for ijk=1:length(ND)
        temptrigs=[temptrigs;ND(ijk).TriggerTimes,...
            ijk*ones(size(ND(ijk).TriggerTimes)),ND(ijk).FFVals,...
            ND(ijk).AmpVals,ND(ijk).RepeatNumber];
    end
    
    if ~isempty(temptrigs)
        [ysort,isort]=sort(temptrigs(:,1));
        temptrigs=temptrigs(isort,:);
    end
    set(handles.TriggerTable,'Data',temptrigs);
end


handles.TAFVALS=handles.TAFVALSALL(handles.NDIndex).vals;
vals=handles.TAFVALS;
guidata(hObject,handles);

v=axis;
axes(handles.ValAxes);hold off;
tmpplt=plot([1:length(vals)].*(2*handles.TemplLen)./fs + handles.TBEFORE,vals,'.-');

ylim(v(3:4));hold on;

handles.ValsPltHndl = tmpplt;
t=handles.t;
tmpplt2=zeros([length(handles.ND(handles.NDIndex).CntRng),1]);
for ii = 1:length(handles.ND(handles.NDIndex).CntRng)
    clr=get(tmpplt(ii),'Color');
    tmpplt2(ii)=plot([t(1),t(end)],[1,1]*handles.ND(handles.NDIndex).CntRng(ii).TH);hold on;
    set(tmpplt2(ii),'Color',get(tmpplt(ii),'Color'));
end
handles.ThreshLineHand = tmpplt2;
handles.TAFVALS=vals;
guidata(hObject,handles);

axes(handles.LabelAxes);

%%%%%%%%%%%%%%%%%
trigtimes = handles.ND(handles.NDIndex).TriggerTimes;
handles.TrigT=trigtimes;
guidata(hObject,handles);
return;

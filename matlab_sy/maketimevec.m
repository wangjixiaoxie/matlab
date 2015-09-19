%output of this function should be a simple 2x1 vector of jesusdates
%7.7.08 updated so as to skip indices with nothing entered.
function [tmvc] = maketimevec(avls)
tmn=avls.tmon;
tmf=avls.tmoff;
tmvc=zeros(length(tmn),2);
for btind=1:length(tmn)
%if date is given, use it.    
%otherwise get the date by finding the corresponding fn.
%Peek at the corresponding index to the batch file.
   if(~isempty(tmn{btind}))
    %open batch file....
        strcmd=['cd '  avls.pvls{btind}];
        eval(strcmd);
        bt=avls.cvl{btind}
        ff=load_batchf(bt);
            fn{1}=ff(1).name;
            disp(fn);
            fn{2}=[];
            indvl=length(ff);
            while (isempty(fn{2}))
                fn{2}=ff(indvl).name;
                indvl=indvl-1;
            end
           

            for ii=1:2
                [pth{ii},nm{ii},ext{ii}]=fileparts(fn{ii});
            
            
            ind{ii}=find(nm{ii}=='_');
            dtvc{ii}=nm{ii}(ind{ii}(1)+1:ind{ii}(2)-1)
            
           
    
            mo(ii)=str2num(dtvc{ii}(3:4));
   
            dy(ii)=str2num(dtvc{ii}(1:2));
            yr(ii)=str2num(['20' dtvc{ii}(5:6)]);
            end
            
            ind_on=find(tmn{btind}==':');
    ind_off=find(tmf{btind}==':');
    %if 7 was entered with no colon
if(isempty(ind_on))   
        hron=str2num(tmn{btind});
        minon=0;
else
    hron=str2num(tmn{btind}(1:ind_on(1)-1));
    minon=str2num(tmn{btind}(ind_on(1)+1:ind_on(1)+2))
end
%if 21was entered with no colon
if(isempty(ind_off))
    hroff=str2num(tmf{btind});
    minoff=0;
else
    hroff=str2num(tmf{btind}(1:ind_off(1)-1));
    minoff=str2num(tmf{btind}(ind_off(1)+1:ind_off(1)+2));
end
if isempty(tmn{btind})|isempty(tmf{btind})
    tmvc=[];
    return;
end
%now turn this date into
%if this is a dirfile
if(~isempty(avls.diron))
    indvl2=find(avls.diron==btind);
    if(indvl2)
        ii=2;
    else
        ii=1;
    end
        
else
    indvl2=[];
    ii=1;
end

tmvc(btind,1)=datenum(yr(1),mo(1),dy(1),hron,minon,0)
tmvc(btind,2)=datenum(yr(ii),mo(ii),dy(ii),hroff,minoff,0)
   
   end
end
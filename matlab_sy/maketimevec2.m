%this function takes a batch file and pulls out the time of the first and
%last file.
function [tmout] = maketimevec2(bt)

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
                tmvc{ii}=nm{ii}(ind{ii}(2)+1:ind{ii}(2)+5)
                mo(ii)=str2num(dtvc{ii}(3:4));
                dy(ii)=str2num(dtvc{ii}(1:2));
                yr(ii)=str2num(['20' dtvc{ii}(5:6)]);
                hr(ii)=str2num(tmvc{ii}(1:2));
                min(ii)=str2num(tmvc{ii}(3:4));
            end
            

tmout(1)=datenum(yr(1),mo(1),dy(1),hr(1),min(1),0)
tmout(2)=datenum(yr(2),mo(2),dy(2),hr(2),min(2),0)
   
  

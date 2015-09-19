function cleandir4(batch,TH,wind,numwind,numnote,CHANSPEC)
% cleandir4(batch,TH,wind,numwind,numnote,CHANSPEC)
% TH is the segment threshold
% wind is the window size in MS
% numwind is the number of notes in a window
% numnote is the numer of times numwind has to be observed in 1 file
%>> cleandir4('batch',1.0e-4,400,4,3);
%>> cleandir4('batch',1.0e-4,400,4,5);
% worked well on bk37w86 :
%>>OLD   ??? cleandir4('batch25',5.0e-5,500,6,10);
%cleandir4('batch',1e4,500,6,10,'obs0');


if (~exist('CHANSPEC'))
    CHANSPEC='obs0r';
end

fid=fopen(batch,'r');
fkeep=fopen([batch,'.keep'],'w');
fdcrd=fopen([batch,'.dcrd'],'w');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	disp(fn);
    goodfile=1;
	
    [pth,nm,ext]=fileparts(fn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn,'0r');
        sm=evsmooth(dat,fs,0.01);
    else
        rd=readrecf(fn);
        if isempty(rd)==0
            [dat,fs]=evsoundin('',fn,CHANSPEC);
                if isempty(dat)==0&(rd.nsamp==length(dat(:,1)))
                    sm=evsmooth(dat,fs,100);
                
                else
                    goodfile=-1;
                end
        else
            goodfile=-1
        end
    end
    keepit=0;
    if(goodfile>0)
        [ons,offs]=evsegment(sm,fs,5.0,30.0,TH);

      
        for ii = 1:length(ons)
            p = find(abs(ons-ons(ii))<=wind);
            if (length(p)>=numwind)
                keepit=keepit+1;
            end
        end
    end
	if(goodfile)
        if (keepit>=numnote)
        
		fprintf(fkeep,'%s\n',fn);
        
        else
		fprintf(fdcrd,'%s\n',fn);
        end
    end
end

fclose(fid);fclose(fkeep);fclose(fdcrd);
return

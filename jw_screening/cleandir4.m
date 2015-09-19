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
% >>OLD? cleandir4('batch',1e4,500,6,10,'obs0');
%
% FOR CATCHING even single FEM CALLS 11/1/12: cleandir4('batch',2.0e-6,400,1,3,'w')
% and for males: cleandir4('cldtest',2.0e-7,500,3,8,'w')
% works 1/2010 (wavread16) cleandir4('batch',5.0e-5,500,6,10,'w')


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
	
    [pth,nm,ext]=fileparts(fn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn,'0r');
        sm=evsmooth(dat,fs,0.01);
    else
        [dat,fs]=evsoundin('',fn,CHANSPEC);
        sm=evsmooth(dat,fs,100);
    end
    [ons,offs]=evsegment(sm,fs,5.0,30.0,TH);
	keepit=0;
	for ii = 1:length(ons)
		p = find(abs(ons-ons(ii))<=wind);
		if (length(p)>=numwind)
			keepit=keepit+1;
		end
	end
	if (keepit>=numnote)
		fprintf(fkeep,'%s\n',fn);
		fprintf(1,'.. was kept. next=')
	else
		fprintf(fdcrd,'%s\n',fn);
		fprintf(1,'.. was dcrd. next=')
	end
end
fclose(fid);fclose(fkeep);fclose(fdcrd);
return

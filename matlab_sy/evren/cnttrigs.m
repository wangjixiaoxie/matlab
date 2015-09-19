function cnt=cnttrigs(batchfile,SHOWTRIGFN,WRTTRIGFNBATCH,USEX);
% cnt=cnttrigs(batchfile,SHOWTRIGFN,WRTTRIGFNBATCH,USEX);
if (~exist('USEX'))
    USEX=0;
else
    if (length(USEX)==0)
        USEX=0;
    end
end

if (~exist('WRTTRIGFNBATCH'))
    WRTTRIGFNBATCH=0;
else
    if (length(WRTTRIGFNBATCH)==0)
        WRTTRIGFNBATCH=0;
    end
end

if (~exist('SHOWTRIGFN'))
    SHOWTRIGFN=0;
else
    if (length(SHOWTRIGFN)==0)
        SHOWTRIGFN=0;
    end
end

cnt=0;
fid=fopen(batchfile,'r');

if (WRTTRIGFNBATCH==1)
    fid2=fopen([batchfile,'.trig'],'w');
    fid3=fopen([batchfile,'.notrig'],'w');
end
if (fid==-1)
	disp(['Batch file : ',batchfile,' does not exist'])';
	return;	
end

while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	rd=readrecf(fn,USEX);
    pp=findstr(fn,'.cbin');
    if (~exist([fn(1:pp(1)),'rec']))
        continue;
    end
	cnt=cnt+length(rd.ttimes);
    if (SHOWTRIGFN)
        if (length(rd.ttimes)>0)
            disp(fn);
        end
    end
    if (WRTTRIGFNBATCH==1)
        if (length(rd.ttimes)>0)
            fprintf(fid2,'%s\n',fn);
        else
            fprintf(fid3,'%s\n',fn);
        end
    end
end
fclose(fid);

if (WRTTRIGFNBATCH==1)
    fclose(fid2);fclose(fid3);
end

disp(['Num of trigs = ',num2str(cnt)]);
return;

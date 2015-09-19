function randsamp(batch,p);
%randsamp(batch,p);

fin=fopen(batch,'r');
fout=fopen([batch,'.rand'],'w');
fout2=fopen([batch,'.notrand'],'w');
while (1)
	fn=fgetl(fin);
	if (~ischar(fn));break;end;
	
	if (~exist(fn,'file')); continue;end;
	
	if (rand(1)<=p)
		fprintf(fout,'%s\n',fn);
    else
        fprintf(fout2,'%s\n',fn);
    end
end
fclose(fin);fclose(fout);fclose(fout2)
return;

function dispwchnt(fv,pp,NT,MKBT); 
%dispwchnt(fv,pp,NT); 
for ii=1:length(pp)
	lbl=fv(pp(ii)).lbl;
	ppp=findstr(lbl,NT);
	pppp=find(fv(pp(ii)).ind==ppp);
	disp([fv(pp(ii)).fn,' ',num2str(fv(pp(ii)).ind),' ',num2str(pppp)]);
end
if (MKBT==1)
	fid=fopen('QBATCHTEMP','w');
	fprintf(fid,'%s\n',fv(pp(1)).fn);
	for ii=2:length(pp)
		if (~strcmp(fv(pp(ii-1)).fn,fv(pp(ii)).fn))
			fprintf(fid,'%s\n',fv(pp(ii)).fn);
		end
	end
	fclose(fid);
end
			
return;

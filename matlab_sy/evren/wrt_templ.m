function wrt_templ(TEMPLNAME,templ);
% wrt_templ(TEMPLNAME,templ);
%
fid=fopen(TEMPLNAME,'w');
for ii=1:size(templ,1)
	fprintf(fid,'%.5e',templ(ii,1));
	for jj=2:size(templ,2)
		fprintf(fid,' %.5e',templ(ii,jj));
	end
	fprintf(fid,'\n');
end
fclose(fid);
return;

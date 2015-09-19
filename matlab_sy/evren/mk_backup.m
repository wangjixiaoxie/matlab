function mk_backup(whichdir,birdname,DOMKDIR,DONOTMAT,DOEBIN);
% mk_backup(whichdir,birdname,DOMKDIR,DONOTMAT,DOEBIN);
%

if (~exist('DONOTMAT','var'))
	DONOTMAT=0;
elseif (length(DONOTMAT)==0)
	DONOTMAT=0;
end
if (~exist('DOMKDIR','var'))
	DOMKDIR=0;
elseif (length(DOMKDIR)==0)
	DOMKDIR=0;
end
if (~exist('DOEBIN','var'))
	DOEBIN=0;
elseif (length(DOEBIN)==0)
	DOEBIN=0;
end

origdir=pwd;disp(origdir);
cd(whichdir);
if (DOEBIN==0)
	if (DONOTMAT==0)
		ff=dir('*.cbin');
	else
		ff=dir('*.cbin.not.mat');
	end
else
	if (DONOTMAT==0)
		ff=dir('*.ebin');
	else
		ff=dir('*.ebin.not.mat');
	end
end

dn=get_dates(ff);

cdir=pwd;
pp=findstr(cdir,birdname);
cdir=cdir((pp(end)+length(birdname)+1):end);

if (~exist([origdir,'/mk_',lower(birdname),'bckups'],'file'))
	INITSTR=1;
else
	INITSTR=0;
end
fid=fopen([origdir,'/mk_',lower(birdname),'bckups'],'a');
if (INITSTR==1)
	fprintf(fid,'#!/bin/sh\n\n');
end
fprintf(fid,'\n%s\n',['cd ',pwd]);
if (DOMKDIR==1)
	str=['mkdir /root/backupdrive/',birdname,'/',cdir,'/'];
	fprintf(fid,'%s\n',str);
end
str=['cp batch* /root/backupdrive/',birdname,'/',cdir,'/'];
fprintf(fid,'%s\n',str);
str=['cp BATCH* /root/backupdrive/',birdname,'/',cdir,'/'];
fprintf(fid,'%s\n',str);
for ii=1:length(dn)
	dv=datevec(dn(ii));
	dv=dv(1:3);
	if (dv(3)<10)
		str=['0',num2str(dv(3))];
	else
		str=num2str(dv(3));
	end
	if (dv(2)<10)
		str=[str,'0',num2str(dv(2))];
	else
		str=[str,num2str(dv(2))];
	end
	tmp=num2str(dv(1));
	str=[str,tmp(3:4)];
		
	str2=['tar -czvf /root/backupdrive/',birdname,'/',...
             cdir,'/',birdname,'_',str,'.tar.gz ',birdname,'_',...
             str,'_*'];
	fprintf(fid,'%s\n',str2);
end
%fprintf(fid,'%s\n',['cd /root']);
%fprintf(fid,'%s\n',['umount backupdrive']);
fclose(fid);
cd(origdir);

function shlabels(batch)

bfile=fopen(batch,'r');
while (~feof(bfile))
	file=fscanf(bfile,'%s',1);
	if file(length(file)-7:length(file))=='.not.mat'
	%if file~=''
		if (~exist(file,'file'))
			disp(['File not found: ',file])
		else
			load(file);
			[r,c]=size(labels);
			if r~=1
			labels=transpose(labels);
			end
			fprintf(1,['File ',file,':\n    ',labels,'\n'])
		end
	%end
	else
		if (~exist([file, '.not.mat'],'file'))
			disp(['File not found: ',file,'.not.mat'])
		else
			load([file,'.not.mat']);
			[r,c]=size(labels);
			if r~=1size(labels)
			labels=transpose(labels);
			end
			fprintf(1,['File ',file,':\n    ',labels,'\n'])
		end
	end
end
fclose(bfile);

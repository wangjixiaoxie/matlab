function lv_randsample_batch(batchkeep,numsongs,starttime)

fid=fopen(batchkeep,'r');
fkeep=fopen([batchkeep '.lvrand.keep'],'w');
fdcrd=fopen([batchkeep '.lvrand.ignore'],'w');

[~,numlines] = fscanf(fid,'%s');
frewind(fid);

for i=1:numlines
    fn{i}=fgetl(fid);
end


time = cellfun(@(x) str2num(x(17:18)),fn);
nlatesongs = sum(time>=starttime);

keepidx = zeros(numlines,1);
if nlatesongs>=numsongs
    disp(['...sampling random songs after ' num2str(starttime)])
    %random sample numsongs
    tempi = randperm(nlatesongs);
    takeidx = tempi(1:numsongs) + sum(time<starttime);
    keepidx(takeidx) =1;
elseif nlatesongs<numsongs & numlines>numsongs
    addfiles = numsongs-nlatesongs;
    disp(['...adding ' num2str(addfiles) ' earlier songs']);
    keepidx(end-nlatesongs-addfiles+1:end)=1;
else
    disp('...keeping all songs')
    keepidx = ones(numlines,1);
end

for i = 1:numlines
    if keepidx(i)
        fprintf(fkeep,'%s\n',fn{i});
    else
        fprintf(fdcrd,'%s\n',fn{i});
    end
end

fclose(fid);
fclose(fkeep);
fclose(fdcrd);

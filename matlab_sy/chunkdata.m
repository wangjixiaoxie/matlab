function [chunkout]=chunkdata(avls,graphvals, offset);
    %start off with time bin as 1 hour
    %1/24 days
    tbin=offset/24
    for ii=1:length(avls.NT)
    for jj=1:length(graphvals.timon)
        
        btind=avls.analind(jj);
        if(~isempty(graphvals.date{btind}));
            vals=getvals(avls.fvcomb{ii},1,'TRIG');
            st=datenum(graphvals.date{btind} ,'yyyy-mm-dd');
            ind=find(floor(vals(:,1))==st);
            chunk=vals(ind,:);
        
        %now make a vector with n timebins and go through finding which
        %indices of the chunk fit into these timebins and calculate the
        %mean and standard deviation of the chunk, output into chunkout
        %
            tminds=[0:tbin:1]
            for kk=1:length(tminds)-1
                tmvl=tminds(kk)
                bastime=floor(chunk(1,1));
                chunkind=find(chunk(:,1)>=bastime+tminds(kk)&chunk(:,1)<bastime+tminds(kk+1))
                if(~isempty(chunkind))
                    chunkout{ii}.mean{jj}(kk)=mean(chunk(chunkind,2));
                    chunkout{ii}.std{jj}(kk)=std(chunk(chunkind,2));
                    chunkout{ii}.tmvl{jj}(kk)=bastime+tminds(kk);
                else
                    chunkout{ii}.mean{jj}(kk)=0
                    chunkout{ii}.std{jj}(kk)=0;
                    chunkout{ii}.tmvl{jj}(kk)=bastime+tminds(kk);
                end
                end
        end
    end
    end
    
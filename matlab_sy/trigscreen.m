function [] = trigscreen(batchfile,thresh,pltboth);
% trigscreen(batchfile,thresh);
%
% goes through the batch file of obs1r files looking for 
% BirdTAF triggers
%

ftrigger = fopen([batchfile,'.trigs'],'a');

if (~exist('pltboth'))
    pltboth=0;
end
if (~exist('thresh'))
	thresh = 300.0;
end

fid = fopen(batchfile,'r');
while (~feof(fid))
	fn  = fscanf(fid,'%s',1);
	disp(fn);
	if (exist(fn,'file'))
		[rs,fs]=evsoundin('',fn,'obs0r');
		[sm,sp,f,t]=evsmooth(rs,fs,100.0);
		sp = log(abs(sp));
        if (pltboth)
            subplot(211);
            title(normstr(fn));
            imagesc(t,f,sp);m=colormap('gray');colormap(m(end:-1:1,:));
            set(gca,'YD','n');axis([0 max(t) 0 8e3]);
            hold on;
        else
            imagesc(t,f,sp);m=colormap('gray');colormap(m(end:-1:1,:));
            set(gca,'YD','n');axis([0 max(t) 0 8e3]);
            title(normstr(fn));
            hold on;
        end
            

		[taf_out,fs]=evsoundin('',fn,'obs1r');
		tt = [1:length(taf_out)]/fs;
		p = find((taf_out(1:end-1)<=thresh)&(taf_out(2:end)>thresh));
        if (pltboth)
            subplot(211);
        end
        hold on;v=axis;
        for ii = 1:length(p)
            plot([1,1]*tt(p(ii)),v(3:4),'r-');
        end
		hold off;
        
        if (pltboth)
            subplot(212);plot(tt,taf_out,'r');
            axis([0 max(t) min(taf_out) max(taf_out)]);
        end
		drawnow;pause;
	end
end
fclose(fid);
return;

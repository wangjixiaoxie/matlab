function [wmeans,wstdv,distrs]=excur(data)
%[means,stdv]=excur(data);

shuff=zeros(length(data),1);
randindex=randperm(length(data));
randindex=randperm(length(data));
for r=1:length(data);
	shuff(r)=data(randindex(r));
end
figure;plot(data,'r');hold on;plot(shuff,'g');

olap=.6;
min_w=5; %may not work with small window lengths (see the for loops). 5 is safe
%max_w=length(data);

wmeans=zeros(length(data),1); %both these vars have zeros for the first few vals (wlength <min_w) oh well it allows the index into these vectors to be meaningful.
wstdv=zeros(length(data),1);
wsterr=zeros(length(data),1);
sdiff=zeros(length(data),1);
distrs=zeros(9,length(data));
%for w=5:100
for w=min_w:length(data)-1;
        count=0;
        incr=floor((1-olap)*w);
        vals=zeros(floor(length(data)/incr),1); %will use count to refer to the indices that get changed only
	valshuff=vals;
        for seg=1:floor(length(data)/incr); %not all of these will make it
                count=count+1;
                %there's at least 2 calcs you can do per window.
                % one is find the range of the data -- if a signal has higher freq content, you get closer to full range per window with smaller window sizes. If it has more LF content, then you need larger window sizes to approach full-range. The other calc to do is abs of the diff b/w initial and final val per window. there will be a peak where window size best "captures" the data's variability.
                if seg>1
		if (incr*seg+w)>length(data);
                        %this may be questionable - but don't want to omit the last few data points. 
                        %first is calc1. 
                        %vals(count)=abs((max(data((length(data)-w):length(data)))-(min(data((length(data)-w):length(data))))));
                        %or choose calc2:
                        vals(count)=abs(data(length(data))-data(length(data)-w));
                        valshuff(count)=abs(shuff(length(shuff))-shuff(length(shuff)-w));
			break;
                end
		end
                %same 2 possible types of calculations here:
                %vals(count)=abs((max(data(incr*seg:(incr*seg+w))))-(min(data(incr*seg:(incr*seg+w)))));
                vals(count)=abs(data(1+incr*(seg-1))-data(incr*(seg-1)+w));
		valshuff(count)=abs(shuff(1+incr*(seg-1))-shuff(incr*(seg-1)+w));
        end
        %get answer for that windowlength here
	distrs(:,w)=hist(vals(1:count),[.5 1 2 4 8 16 32 64 128]); %to preserve full distributions per windowlength
	%or, just extract p-values for identity of shuff and nonshuff:
	[ht(w),pt(w)]= ttest2(vals(1:count),valshuff(1:count),'Vartype','unequal','Tail','right'); % right tail means test equality against alt hyp that first vector's pop mean is > second vector's pop mean

	[hks(w),pks(w)] = kstest2(vals(1:count),valshuff(1:count),'Tail','smaller'); % smaller means test equality against alt hyp that data vals are larger than shuff vals ("smaller" because cdf of data is smaller than for shuff if data vals tend larger than shuff's / rightward shift)

	if w<300
	if w>296
%		distrs(:,w)		
%		count
%		figure;bar(distrs(:,w));title([w,' data'])
%		distrshuff(:,w)=hist(valshuff(1:count),[.5 1 2 4 8 16 32 64 128]);figure;bar(distrshuff(:,w));title([w,' shuff']);
	end
	end
	wmeans(w)=mean(vals(1:count));
        wstdv(w)=std(vals(1:count));
	wsterr(w)=wstdv(w)./(sqrt(count));
	wshuffmeans(w)=mean(valshuff(1:count));
	sdiff(w)=wmeans(w)-wshuffmeans(w);
        %delete this one - just for checking that counts per wlength are reasonable:
        %counts(w)=count;
end
figure;
plot(pt);title(inputname(1));%hold on;plot(pks,'r');
%figure;
%figure;plot(1:1000,wmeans,'o',1:1000,wmeans+wsterr,1:1000,wmeans-wsterr);

%current 2 lines output plots:
%errorbar(wmeans,wsterr);hold on;plot(wmeans,'r');plot(sdiff,'g');
%surf(distrs);

%plot(wmeans+wsterr);hold on;plot(wmeans-wsterr);
%plot(wmeans+wstdv);hold on;plot(wmeans-wstdv);

%figure;plot(counts); %delete this line if is reliably a smooth decay
return

function interval=confid(dat,intv);
%interval=confid(dat,intv);

dat2=sort(dat);
mm=median(dat);
[y,i]=min(abs(dat2-mm));
ind=[i,i];

nn=1;ln=length(dat);
mn=dat2(ind);mx=mn;
while (1)
	if (ind(1)<1)|(ind(2)>length(dat2))
		break;
	end

	ind = ind+[-1,1];
    if (ind(1)<1)
        ind(1)=1
    end
    if (ind(2)>length(dat2))
        ind(2)=length(dat2)
    end
	mn=dat2(ind(1));mx=dat2(ind(2));
	pp=find((dat2>=mn)&(dat2<=mx));

	if (length(pp)>=intv*ln)
		break;
    end
    if (ind(1)==1)&(ind(2)==length(dat2))
        break;
    end
end
interval=[mn,mx];
return;

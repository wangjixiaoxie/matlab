figure
recbs;
sumdata=[];
for ii=1:length(recbs)
    LW=1;
    rb=recbs(ii)
    shmag=rb.sh-rb.bas;
    if(shmag<0)
        drxn=2;
    else
        drxn=1;
    end
    tst(1)=1;;
    if(drxn==1)
        tst(2:4)=(rb.rec-rb.bas)./abs(shmag);
    else
        tst(2:4)=(rb.bas-rb.rec)./abs(shmag);
    end
    if(~isempty(rb.mu))
        outvl=abs(rb.mu-rb.bas)/shmag;
        plot(1,outvl,'ro');
        hold on;
        if(ii<3)
        plot([1.2 1.2],[(outvl-rb.std/shmag) (outvl+rb.std/shmag)])
        else
              plot([1 1],[(outvl-rb.std/shmag) (outvl+rb.std/shmag)])
        end
              LW=2;
    end

    plot(1:4,tst,'ko');
    hold on;
    plot(1:4,tst,'k','Linewidth',LW)
    hold on;
    sumdata=[sumdata;tst]
end
    
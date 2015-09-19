function deltas=goodrep
window=[32,64,128];
shiftover=[4,8,16];
target=jc_flucRW;
target=target/(max(abs(target))+1);
signal=jc_fluca(target);
target=target*100+2300;
deltas=zeros(3,3);
for i=1:3
    windowsize=window(i);
    for j=1:3
        shift=shiftover(j);
        st=jc_evstftoneharmonic(signal,2000,2600,windowsize,shift);
        resampled_target=resample(target,1,shift);
        kk=1;
        shifter=(windowsize/shift)/2;
        del=0;
        for ii=20:180
            del=del+(st(ii)-resampled_target(ii+shifter)).^2;
            kk=kk+1;
        end
        deltas(i,j)=del/kk;
    end
end

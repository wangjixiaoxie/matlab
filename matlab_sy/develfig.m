%script for analyzing juvenile song


datin{1}=rawsong(floor(3.6*fs):floor(5.2*fs));
datin{2}=rawsong3(floor(2.8*fs):floor(4.4*fs));
datin{3}=rawsong2(floor(1.3*fs):floor(2.9*fs));


for ii=1:3
subplot (3,1,ii)

 [sm,sp,t,f]=evsmooth(datin{ii},fs,.005);
 imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
end

rawong=629.104

rawsong2=711.66
rawsong3=705.131

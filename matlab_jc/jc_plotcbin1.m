figure
clear fnvl
clear bnds
%pathvl{1}='/oriole/pk32bk28/ac325'
fnvl{1}='50_395527132_4_14_17_7_2.wav'
bnds{1}=[10 13.5]
% pathvl{2}=pathvl{1}
% fnvl{2}='r95pk42_270507_1051.8688.cbin'
% bnds{2}=[3 6]
% 
% 
% pathvl{3}='/doyale/twarren/r95pk42/templteste3/'
% fnvl{3}='o62bk75_210507_0911.2654.cbin'
% bnds{3}=[2.7 5.7]
% pathvl{4}=pathvl{3}
% fnvl{4}='o62bk75_210507_1216.4231.cbin'
% bnds{4}=[2.5 5.5]

figure
for ii=1:length(fnvl)
subplot(length(fnvl),1,ii)
strcmd=['cd ' pathvl{ii}]
eval(strcmd);
[plainsong,fs]=wavread('50_395527132_4_14_17_7_2.wav');
%here for wav, =wavread('');
plainsong=floor(plainsong(bnds{ii}(1)*fs:bnds{ii}(2)*fs));
[sm,sp,t,f]=evsmooth(plainsong,fs,30);
 imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
 
 colormap('hot')
 box off;
 axis on;   

end
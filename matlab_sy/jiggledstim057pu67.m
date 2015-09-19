%comments on how this algorithm works%
%The mask is a transition, balancing 1/2 how the first bit of data is used, and 
% 1/2 how the second bit of data is used
%



%This will be used to make random pitch shifts

%First one pitch shift
    
%37 syllables in song.  Divide into 12,12,13

%randvec=floor(3*rand(length(xlist)/2,1));
randvec(1:13)=0;
randvec(14:25)=1;
randvec(26:37)=2;
permout=randperm(37)
randvec=randvec(permout);
    maskwid=.005
    
    datnorm_use=corpshiftednormg{1};
    datup_use=corpshiftednormg{3};
    datdown_use=corpshiftednormg{2};
    dat_shift=datnorm_use;


for i=1:2:length(xlist)
    if(randvec(round(i/2))==0)
        switchsong=datnorm_use;
    elseif(randvec(round(i/2))==1)
        switchsong=datdown_use;
    else
        switchsong=datup_use;
    end
        
    
    
    tlim2=[xlist(i),xlist(i+1)];
    
    mask=mk_mask(datnorm_use,fs,tlim2(1),tlim2(2),maskwid);
    dat_shift=dat_shift.*(1.0-mask) + switchsong.*mask;
    %[sm,sp,t,f]=evsmooth(dat_upsngl,fs,0.01);
    %[smnormu,spnormu,tnormu,fnormu]=evsmooth(datnorm_use,fs,0.01);
    %figure(3);clf;ax=zeros([2,1]);
    %ax(1)=subplot(211);imagesc(tnormu,fnormu,log(abs(spnormu)));syn;ylim([0,1e4]);
    %grid on;
    %ax(2)=subplot2(212);imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    %grid on;
    %linkaxes(ax);
    
end

%To plot this,need to create a vector sampled every 5ms.
%Start with a vector of zeros
%for each entry in xlist, find the correct points and change them.
clear imagevec
totaltime=length(dat_shift)/fs;
imagesr=1000;
imagevec=0:1/imagesr:totaltime;
imaget=imagevec;
imagevec(:)=-.2;



for i=1:2:length(xlist)
    imagevec(floor(xlist(i)*imagesr):floor(xlist(i+1)*imagesr))=randvec(round(i/2));
end

figure;
ax(1)=subplot(10,2,2:2:16)
[sm,sp,t,f]=evsmooth(dat_shift,fs,0.01);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
%plot(0:1/fs:(totaltime-1/fs),dat_shift)
ax(2)=subplot(10,2,18:2:20)
imagesc(imaget,f,imagevec)
ax(3)=subplot(10,2,1:2:15)

%This is the end of what I used.



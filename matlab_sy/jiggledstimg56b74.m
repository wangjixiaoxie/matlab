%comments on how this algorithm works%
%The mask is a transition, balancing 1/2 how the first bit of data is used, and 
% 1/2 how the second bit of data is used
%



%This will be used to make random pitch shifts

%First one pitch shift
    
%37 syllables in song.  Divide into 12,12,13

%randvec=floor(3*rand(length(xlist)/2,1));
%measure syllables in song
figure;
subplot(1,1,1)
   ax(1)=subplot(2,1,1)
    [sm,sp,t,f]=evsmooth(corpshiftednormg{1},fs,0.01);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
   ax(2)=subplot(2,1,2)
        plot(0:1/fs:(length(corpshiftednormg{1})-1)/fs,corpshiftednormg{1})
   
linkaxes(ax,'x');
[x,y]=ginput();


randvec(1:16)=0;
randvec(17:33)=1;

permout=randperm(33)
randvec=randvec(permout);
    maskwid=.005
    
    datnorm_use=corpshiftednormg{1};
    datup_use=corpshiftednormg{3};
    datdown_use=corpshiftednormg{2};
    dat_shift=datnorm_use;


  %2nd iteration  
for i=1:length(randvec)
    if(randvec(i)==0)
        switchsong{1}=datnorm_use;
        switchsong{2}=datdown_use;
    elseif(randvec(i)==1)
        switchsong{1}=datdown_use;
        switchsong{2}=datnorm_use;
    else
        switchsong=datup_use;
    end
        
    
    
    tlim2=[xlist(i),xlist(i+1)];
    
    mask=mk_mask(datnorm_use,fs,tlim2(1),tlim2(2),maskwid);
    
    %dat_shift is the mask, to keep previous pitch changes
    dat_shift=dat_shift.*(1.0-mask) + switchsong{1}.*mask;
    dat_shift2=dat_shift2.*(1.0-mask) + switchsong{2}.*mask;
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
imagevec(:)=-1;



for i=1:length(xlist)-1
    imagevec(floor(xlist(i)*imagesr):floor(xlist(i+1)*imagesr))=randvec(i);
end

figure;
ax3(1)=subplot(10,2,2:2:16)
[sm,sp,t,f]=evsmooth(dat_shift,fs,0.01);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
%plot(0:1/fs:(totaltime-1/fs),dat_shift)
ax3(2)=subplot(10,2,18:2:20)
imagesc(imaget,f,imagevec)
ax3(3)=subplot(10,2,1:2:15)
[sm,sp,t,f]=evsmooth(dat_shift2,fs,0.01);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
linkaxes(ax3)
%This is the end of what I used.



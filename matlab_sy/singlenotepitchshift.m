%comments on how this algorithm works%
%The mask is a transition, balancing 1/2 how the first bit of data is used, and 
% 1/2 how the second bit of data is used
%



%This will be used to make random pitch shifts

%First one pitch shift
    
%37 syllables in song.  Divide into 12,12,13

%randvec=floor(3*rand(length(xlist)/2,1));

tlim2(1)=xlist(1);
tlim2(2)=xlist(2);


    maskwid=.005
    
    datnorm_use=corpshiftednormg{1};
    datup_use=corpshiftednormg{3};
    datdown_use=corpshiftednormg{2};
    

    mask=mk_mask(datnorm_use,fs,tlim2(1),tlim2(2),maskwid);
    datup_single=datnorm_use.*(1.0-mask) + datup_use.*mask;
    datdown_single=datnorm_use.*(1.0-mask) + datdown_use.*mask;
    %[sm,sp,t,f]=evsmooth(datup_single,fs,0.01);
    %[smnormu,spnormu,tnormu,fnormu]=evsmooth(datdown_single,fs,0.01);
    %figure(3);clf;ax=zeros([2,1]);
    %ax(1)=subplot(211);imagesc(tnormu,fnormu,log(abs(spnormu)));syn;ylim([0,1e4]);
    %grid on;
    %ax(2)=subplot(212);imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    %grid on;
    %linkaxes(ax);
    


%To plot this,need to create a vector sampled every 5ms.
%Start with a vector of zeros
%for each entry in xlist, find the correct points and change them.
wavwrite(datup_single,fs,'5767single_up.wav');
wavwrite(datdown_single,fs,'5767single_down.wav');
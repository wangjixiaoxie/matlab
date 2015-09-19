%plotpanel2.m
%this is called by plotstimexample.m
%9.13.09 for stim data

%this also plots the entire time course.

%this is bk48w74
stsind=1
PLOTALL=1
PLOTINIT=0
figure
ax(1)=subplot(121)
ax(2)=subplot(122)

clear ind
cmd=['cd ' sts(stsind).path];
eval(cmd);

cmd=['load ' sts(stsind).matfilename]
eval(cmd);

fb=avls.fbmean
ct=avls.ctmean
sfact=sumbs(stsind).sfact(1);
 ct=ct./sfact;
    fb=fb./sfact;
ind{1}=[1:7 13 ]
ind{2}=[18 19 24 28 29 31 33]
ind{3}=[1:7 13 14 15 17 18]
ind{4}=[18 19 24 25 28 29 31 33 35 36 39 41 42 43]
%fill timevecs

if(PLOTINIT)
    for ii=1:length(ax)
        axes(ax(ii));
        crind=ind{ii};
        rawtms=sumbs(stsind).flrtmvec(crind)  
        plot_tmcourse_general(rawtms,ct(crind),fb(crind));
        axis([-1 9 6700/sfact 7800/sfact])
        axis square
        plot([-1 9],[sumbs(stsind).initmean/sfact sumbs(stsind).initmean/sfact],'k--')
    end
end

if(PLOTALL)
    for ii=1:length(ax)
        axes(ax(ii));
        crind=ind{ii+2}
        rawtms=sumbs(stsind).tmvec(crind);
        plot_tmcourse_general(rawtms,ct(crind),fb(crind));
        axis([-1 14 2.2 2.6])
        axis square
        plot([-1 9],[sumbs(stsind).initmean/sfact sumbs(stsind).initmean/sfact],'k--')
    end
end
        
 




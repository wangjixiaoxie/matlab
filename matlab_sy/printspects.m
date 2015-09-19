
%instead of the files, use the indices to the fvstructure, that way 
%I can print out six seconds centered around the note.




% g21g81_240107_0902.192.cbin
% 185
% 
% g21g81_240107_1910.933.cbin
% 190
% 
% g21g81_260107_1013.193.cbin
% 196
% 
% 
% 2/3
% 6
% 
% 2/10
% 82
% 
% 2/12
% 98
% 103
% 
% 1/13
% 114
% 
% 2/14
% 167
% 
% 2/15
% 172
% 
% 2/19
% 1506

loc{1}='/doyale2/g21g81/postscreen/'
loc{2}='/doyale2/g21g81/10microstim/'
loc{3}='/doyale1/twarren/g21g81/stimoff/'


fv_vls=[185 190 196 6 82 98 103 114 167 172 1506 2704]
locvals=[1 1 1 2 2 2 2 2 2 2 3 3 ]

pltsperpg=4
colormap('hot')

for fvind=1:length(fv_vls)
    curfv=fv_vls(fvind);
    %sets the new figure on 1, 1+pltsperpg....
    if(mod(fvind, pltsperpg)==0)
        figure
    end

    %here 1, 1 +plts perpg should be 1.
    pltindx= mod(fvind, pltsperpg)+1
    subplot(pltsperpg,1,pltindx)
    
    %get in proper directory
    strcmd=['cd ' loc{locvals(fvind)}]
    eval(strcmd);
    
    %now, 
    
    
    [datfile,fs]=evsoundin('',fvcomb(curfv).fn,'obs0');
    %to find time
    
    ontm=fvcomb(curfv).ons(fvcomb(curfv).ind);
    
    start_time=(ontm/1000)*fs -2*fs;
    if (start_time<0)
        start_time=0;
    end
    
    end_time=(ontm/1000)*fs+2*fs;
    if (end_time)>length(datfile);
        end_time=length(datfile);
    end
    
    
    datfile=datfile(start_time:end_time);
    
    [sm,sp,t,f]=evsmooth(datfile,fs,.005);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    box off;
    title(fvcomb(curfv).fn)
    %make this a function, like print spect
    
end
    


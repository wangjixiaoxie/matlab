function cleandir_spect(bname,TH,wind,numwind,numnote,CHANSPEC)
%
% this version (sam sober, 12/13/2005) uses TWO criteria to include or
% reject files.  one is the threshold crossing criterion used in
% the standard cleandir.m programs.  the new criterion i comparing
% the ration of the the mean PSD in the high (>10kHz) and low (<10kHz)
% frequency band.  empirically, the ration of hi:lo is much higher
% (by about 2 orders of magnitude) for files containing only movement
% artifact comparted to files containing only song.
%
%
% % for CBS (wav) files:
% cleandir_spect('bnamefoo',10^-5,500,6,1)
%
% % for cbin files:
% cleandir_spect('bnamefoo',10^6,500,6,1)
%
% name of bname file - names of files  - takes either WAV or .cbin (assumes 'obs0'
% threshold of smoothed waveform (evsmooth - check)
% window size (see below)
% number of crossings (below)
% numnote (below)
%
%  This file produced two output files - bname.keep and bname.dcrd, both
%  containing names of files to be kept or discarded.
%
%  If files in bname are .cbin, then the bname.dcrd will contain the names
%  of all appropriate .rec and .tmp files too.
%
%   NOTE:   this file needs to call evread_obsdata.m, evsmooth.m,
%   evsegment.m, spect_from_waveform.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
%  New feature -regular cleandir_spect.m is slowed down by evsmooth.m,
%  which is used to determine number of thresh crossings and takes about
%  10x longer to run than spect_from_waveform.m.
%
%  The new feature will first check for spectral criterion, and if file
%  fails this, will skip the test for number of crossings.
%1

% if true, displays time spent on various computations
display_time=0;

if (~exist('CHANSPEC'))
    CHANSPEC='obs0';
    %    CHANSPEC='obs0r';
end
save_or_del_code_save=[];
hi_over_lo_save=[];
filename_save=[];
ct=1;
%eval(sprintf('fid=fopen(%s,"r")',bname))
fid=fopen(bname,'r');
fkeep=fopen([bname,'.keep'],'w');
fdcrd=fopen([bname,'.dcrd'],'w');
while (1)
    fn=fgetl(fid);
    %    fn_tmp=[fn(1:end-5) '.tmp'];
    %    fn_rec=[fn(1:end-5) '.rec'];
    if (~ischar(fn))
        break;
    end

    if fn(end-3)=='.'   % .mat, .jpg etc
        fn=fgetl(fid);
    end
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    disp(fn);


    [pth,nm,ext]=fileparts(fn);
    if ~strcmp(fn(end-3:end),'cbin')
        tic
        [dat,fs,NBITS]=wavread(fn);
        tmcnt=toc;
        if display_time
        disp(['Time to load WAV file - ' num2str(tmcnt) ' sec'])
        end
        is_cbin=0;
    else
        tic
        [dat,fs]=evsoundin('',fn,'obs0');
        tmcnt=toc;
        if display_time
        disp(['Time to load cbin file - ' num2str(tmcnt) ' sec'])
        end
       % smoothing is the long step - do this later
%        sm=evsmooth(dat,fs,100);
        is_cbin=1;
        fn_tmp=[fn(1:end-5) '.tmp'];
        fn_rec=[fn(1:end-5) '.rec'];
        fn_neuralnot=[fn(1:end-5) '.cbin.neuralnot.mat'];
    end
    tic
    [S1,F1,T1,P1] =spect_from_waveform(dat,fs,0,[0 16]); % will give a 512-point transform,
    tmcnt=toc;
    if display_time
    disp(['Time to do spect - ' num2str(tmcnt) ' sec'])
    end
    hi_F_id=find(F1>10000);
    lo_F_id=find(F1<=10000);
    mean_pow_hi=mean(mean(P1(hi_F_id,:)));
    mean_pow_lo=mean(mean(P1(lo_F_id,:)));
    hi_over_lo=mean_pow_hi/mean_pow_lo;
    hi_over_lo_save(ct)=hi_over_lo;
    filename_save{ct}=fn;

    
    % if hi_over_lo>.05;  % HIGH ratio indicates hopping around
    if hi_over_lo>.1;  % HIGH ratio indicates hopping around
        keepit_ratio=0;
    else
        keepit_ratio=1;
    end
    
    % THE NEW STUFF - if spectral requirement is met, THEN look at number
    % of crossings.  if not, just set keepit_crossing=0
    if keepit_ratio
        tic
        if ~strcmp(fn(end-3:end),'cbin')
            sm=evsmooth(dat,fs,0.01);
        else
            sm=evsmooth(dat,fs,100);
        end
        tmcnt=toc;
        if display_time
        disp(['Time to smooth file - ' num2str(tmcnt) ' sec'])
        end
        %            figure(100);clf;subplot(2,1,1);plot(log10(sm));subplot(2,1,2);plot(sm);return
        [ons,offs]=evsegment(sm,fs,5.0,10.0,TH);

        keepit_crossing=0;
        for ii = 1:length(ons)
            p = find(abs(ons-ons(ii))<=wind);
            if (length(p)>=numwind)
                keepit_crossing=keepit_crossing+1;
            end
        end
        crossing_shortcut=0;
    else
        crossing_shortcut=1;
        keepit_crossing=0;
    end

    if (keepit_crossing>=numnote) & keepit_ratio
        disp('KEEPING')
        save_or_del_code_save(ct)=0;
        fprintf(fkeep,'%s\n',fn);
        if is_cbin % write tmp and rec to keep
            %fprintf(fkeep,'%s\n',fn_tmp);
            %fprintf(fkeep,'%s\n',fn_rec);
            %            fprintf(fkeep,'%s\n',fn_neuralnot);
        end
    else
        if ~keepit_ratio & crossing_shortcut
            disp('Rejecting on basis of spectral ratio - NOT CHECKING CROSSINGS')
            save_or_del_code_save(ct)=1;
        elseif ~(keepit_crossing>=numnote) & keepit_ratio & ~crossing_shortcut
            disp('Rejecting on basis of number of crossings (SPECT WAS OK)')
            save_or_del_code_save(ct)=2;
        else
            error('hjere')
            disp('Rejecting on basis of both criteria')
            save_or_del_code_save(ct)=3;
        end
        fprintf(fdcrd,'%s\n',fn);
        if is_cbin
            fprintf(fdcrd,'%s\n',fn_tmp);
            fprintf(fdcrd,'%s\n',fn_rec);
            %            fprintf(fdcrd,'%s\n',fn_neuralnot);

        end
    end
    ct=ct+1;
end
save ratio hi_over_lo_save filename_save save_or_del_code_save
fclose(fid);fclose(fkeep);fclose(fdcrd);
return

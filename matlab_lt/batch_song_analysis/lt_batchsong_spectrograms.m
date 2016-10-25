clear all; close all;
songfn='rd12pu6_PBS_170915_135530.8070.cbin';

% =============== Read song file
% 1) read cbin file
[dat, Fs, DOFILT, ext]=ReadDataFile(songfn,'0');

if (1)
    lt_plot_spectrogram(dat, Fs, 1, 0);
else
% 2) smooth, filter sound file
[sm,sp,t,f]=SmoothData(dat,Fs,DOFILT,'‘hanningfirff’'); % ends up doing buttor, with filtfilt. (because there is extra quote marks in filter type)
% default: 512 NFFT, using 512 sample windows, 0.8 olap
% sm: raw data --> filtered --> squared --> smoothed (uniform 2ms window), 

% 3) Since song dat was bandpassed (500 to 10000), remove freqs that don't
% want
sp=sp(10:140,:);
f=f(10:140);

% 4) plot
figure;
% imagesc(t, f, log(sp));
% axis([t(1) t(end) f(1) f(end)]);

% 4) plot
%         pp=find(sp>0);
%         mntmp = min(min(sp(pp)));
%         pp=find(sp==0);
%         sp(pp) = mntmp;
%         
%         % second, take log
%         sptemp=log(sp);
%         sptemp = sptemp - min(min(sptemp));
%         sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
%         
%         
%         % PLOT
%         % Prepare subplot slot
%                 sp=sp(10:140,:);
%         f=f(10:140);
%         
%         sp=sp(:,312:end); % remove 1st sec.
%         t=t(312:end);
%         
        
        
% 4) THIS GOOD - 
        % TAKE LOG of sp
        % first, convert any sp values of 0 to non-zero(to the lowest value present);
        % solves problem of taking log of 0
        pp=find(sp>0);
        mntmp = min(min(sp(pp)));
        pp=find(sp==0);
        sp(pp) = mntmp;
        
        % second, take log
        sptemp=log(sp);
        sptemp = sptemp - min(min(sptemp));
        sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
        
                
        % Plot
        imagesc(t, f, flipud(sptemp));
        
        axis([t(1) t(end) f(end) f(1)]);
        axis([t(1) t(end) f(1) f(end)]);
        end
        
%% ====== PLAY AUDIO

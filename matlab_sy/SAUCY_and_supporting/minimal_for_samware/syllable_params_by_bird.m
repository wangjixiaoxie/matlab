function [f_cutoff,t_assay,spect_params,syl_name_in_notmat]=syllable_params_by_bird(bname,syl)

if strcmp(bname,'W15-W94')% MEL
    if strcmp(syl,'g')
        f_cutoff=[2400 3400];
        t_assay=.016;
        spect_params=[0 8];
    end
elseif strcmp(bname,'ZFg59o19')
    if strcmp(syl,'a') %syl a -
        f_cutoff=[1000 2500];
        t_assay=.008;
        spect_params=[0 16];
    elseif strcmp(syl,'b') % syl b - not great
        f_cutoff=[3500 5500];
        t_assay=.016;
        spect_params=[.5 16];
    elseif strcmp(syl,'c')  % syl c - good for f1, but this is v. low freq.
        %spect. complex at higher freq
        f_cutoff=[1400 1800];
        t_assay=.098;
        spect_params=[0 64];
    elseif strcmp(syl,'d') % syl d - good, though dense stack.  this is the FIFTH harmonic
        %        f_cutoff=[2600 3400];
        f_cutoff=[2800 3400];
        t_assay=.016;
        spect_params=[0 64];
    end

elseif strcmp(bname,'ZFg13w11')
    if strcmp(syl,'c') % old syl b -
        f_cutoff=[1000 2000];
        t_assay=.04;
        spect_params=[0 32];
    elseif strcmp(syl,'d') % old syl c -
        f_cutoff=[2500 3000];
        t_assay=.016;
        spect_params=[0 32];
    elseif strcmp(syl,'b') % old syl c -
        f_cutoff=[1400 2500];
        t_assay=.004;
        spect_params=[0 8];% maybe try [0 16]?
    end
    %     if strcmp(syl,'b') %syl b -
    %         f_cutoff=[1400 1800];
    %         t_assay=.056;
    %         spect_params=[0.8 64];
    %     elseif strcmp(syl,'c') %syl c -
    %         f_cutoff=[1900 2500];
    %         t_assay=.016;
    %         spect_params=[0 32];
    %     end
elseif strcmp(bname,'ZFg5o78')  % trying 5,000,000 thresh, min interval 2
    if strcmp(syl,'a')
        f_cutoff=[2000 6000];
        t_assay=.008;
        spect_params=[0 16];
    elseif strcmp(syl,'b')
        f_cutoff=[2500 4000];
        t_assay=80;% means t_pct=.8
        spect_params=[0 16];
    elseif strcmp(syl,'c')
        f_cutoff=[2800 3600];
        t_assay=.1;
        spect_params=[0 64];
    elseif strcmp(syl,'d')
        f_cutoff=[2600 3400];
        t_assay=0.075;
        spect_params=[0.5 32];
    end

    % elseif strcmp(bname,'ZFg5o78')  % old
    %     if strcmp(syl,'a')
    %         f_cutoff=[2000 6000];
    %         t_assay=.02;
    %         spect_params=[0 8];
    %     elseif strcmp(syl,'b')
    %         f_cutoff=[2500 4000];
    %         t_assay=80;% means t_pct=.8
    %         spect_params=[0 16];
    %     elseif strcmp(syl,'c')
    %         f_cutoff=[2750 3250];
    %         t_assay=.1;
    %         spect_params=[0 64];
    %     elseif strcmp(syl,'d')
    %         f_cutoff=[2750 3250];
    %         t_assay=0.075;
    %         spect_params=[0.5 32];
    %     end

elseif strcmp(bname,'g18g8')
    if strcmp(syl,'a')
        f_cutoff=[600 2750];
        t_assay=40;
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[2000 4000];
        t_assay=.016;
        spect_params=[0.5 8];
        %         f_cutoff=[1000 3000];
        %         t_assay=30;
        %         spect_params=[0 8];
    elseif strcmp(syl,'c')
        f_cutoff=[4000 6000];
        t_assay=.016;
        spect_params=[0.5 16];
    elseif strcmp(syl,'d')
        f_cutoff=[1500 3000];
        t_assay=.04;
        spect_params=[0.5 16];
    elseif strcmp(syl,'e')
        f_cutoff=[3000 5000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'f')
        f_cutoff=[1000 3000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'g')
        f_cutoff=[900 1750];
        t_assay=75;
        spect_params=[0 8];
    elseif strcmp(syl,'i')
        f_cutoff=[2500 4500];
        t_assay=50;
        spect_params=[0 16];
    elseif strcmp(syl,'j')
        f_cutoff=[500 1100];
        t_assay=25;
        spect_params=[0.5 8];
    end
elseif strcmp(bname,'pu24w39')
    if strcmp(syl,'a')
        f_cutoff=[1500 3500];
        t_assay=20; % means t_pct=.2
        spect_params=[0 8];
    elseif strcmp(syl,'e')
        f_cutoff=[1000 3000];
        t_assay=20; % means t_pct=.2
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[3000 5000];
        t_assay=15; % means t_pct=.15
        spect_params=[0 8];
    elseif strcmp(syl,'B')
        f_cutoff=[3000 4750];
        t_assay=.004;
        spect_params=[0 8];
    elseif strcmp(syl,'i')
        f_cutoff=[3000 5000];
        t_assay=.024;
        spect_params=[0 16];
    elseif strcmp(syl,'c')
        f_cutoff=[3000 5000];
        t_assay=.004;
        spect_params=[0 8];
    elseif strcmp(syl,'d')
        f_cutoff=[4000 8000];
        t_assay=.004;
        spect_params=[0 8];
    end
elseif strcmp(bname,'pu44w52')
    if strcmp(syl,'d')
        f_cutoff=[6000 8000];
        t_assay=40;
        spect_params=[0.5 16];
    elseif strcmp(syl,'c') % no pitch corrs with any, the low-pitch params
        % get better corrs with amplitude, but this probably has more to do
        % with use of a pct. - could try looking at higher pitch with a
        % pct. time marker...
        f_cutoff=[3000 6000];
        t_assay=50;
        spect_params=[.5 16];

        % % e.g.could try this one...
        %         f_cutoff=[6000 10000];
        %         t_assay=30;
        %         spect_params=[.5 8];

        % old one used for higher -pitched assay of 'c'
        %         f_cutoff=[6000 10000];
        %         t_assay=.012;
        %         spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[3000 7000];
        t_assay=.008;
        spect_params=[.5 16];
        %        t_assay=.012;
        %       spect_params=[.75 16];
    elseif strcmp(syl,'a')
        f_cutoff=[2700 5000];
        t_assay=25;
        spect_params=[0 8];
    elseif strcmp(syl,'i')
        f_cutoff=[2500 4500];
        t_assay=50;
        spect_params=[0 16];
    elseif strcmp(syl,'j')
        f_cutoff=[500 1100];
        t_assay=25;
        spect_params=[0.5 8];
    end

    % MEL's bird
elseif strcmp(bname,'G26-G23-05282006') | strcmp(bname,'G26-G23-05302006') | strcmp(bname,'G26-G23-05292006')| strcmp(bname,'g26g23')

    if strcmp(syl,'a')
        f_cutoff=[1000 3000];
        %        f_cutoff=[3000 5000];
        t_assay=.004;
        spect_params=[0 8];
    elseif strcmp(syl,'b')  % CANT QUANT PITCH - THIS IS FOR PREMOTOR WINDOW
        f_cutoff=[1000 8000];
        t_assay=.004;
        spect_params=[0 8];
    elseif strcmp(syl,'c')
        %        f_cutoff=[3000 5000];
        f_cutoff=[3500 5000];
        t_assay=.02;
        spect_params=[0 8];
    elseif strcmp(syl,'d')
        f_cutoff=[3500 5000];
        t_assay=.024;
        spect_params=[0.5 16];
    elseif strcmp(syl,'e')
        %        f_cutoff=[3000 5000];
        f_cutoff=[3000 4400];
        t_assay=.004;
        spect_params=[0 8];
    elseif strcmp(syl,'f')
        f_cutoff=[1500 2500];
        t_assay=50;
        spect_params=[0 8];
    elseif strcmp(syl,'g')
        f_cutoff=[6000 8000];
        t_assay=.024;
        spect_params=[0.5 16];
    end
elseif strcmp(bname,'r75g59')   % Jon's bird
    if strcmp(syl,'a')
        f_cutoff=[3000 5000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[3000 6000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'c')
        f_cutoff=[750 2000];
        t_assay=70;
        spect_params=[0 8];
    elseif strcmp(syl,'d')
        f_cutoff=[1000 3000];
        t_assay=70;
        spect_params=[0 8];
    elseif strcmp(syl,'e')
        f_cutoff=[3500 5500];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'f')
        f_cutoff=[4500 6500];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'1')
        f_cutoff=[4500 6500];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'g')
        f_cutoff=[1000 3000];
        t_assay=50;
        spect_params=[0 8];
    elseif strcmp(syl,'h')
        f_cutoff=[500 2000];
        t_assay=95;
        spect_params=[0 8];
    end
elseif strcmp(bname,'pu5b3')   % Jon's bird
    if strcmp(syl,'a')
        f_cutoff=[2000 2500];
        t_assay=75;
        spect_params=[0.5 16];
    elseif strcmp(syl,'b')
        f_cutoff=[3000 5000];
        t_assay=.024;
        spect_params=[0 16];
    elseif strcmp(syl,'g')
        f_cutoff=[3500 5000];
        t_assay=.016;
        spect_params=[0.5 16];
    elseif strcmp(syl,'c')
        f_cutoff=[3000 4600];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'d')
        f_cutoff=[1000 3000];
        t_assay=40;
        spect_params=[0 8];
    elseif strcmp(syl,'e')
        f_cutoff=[1500 3500];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'f')
        f_cutoff=[1000 3500];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'i')
        f_cutoff=[3000 4500];
        t_assay=.024;
        spect_params=[0 16];
    elseif strcmp(syl,'j')
        f_cutoff=[500 1500];
        t_assay=50;
        spect_params=[0 16];
    end
elseif strcmp(bname,'w48o28')   % Jon's bird
    if strcmp(syl,'a')
        f_cutoff=[3000 5000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[1500 2500];
        t_assay=.018;
        spect_params=[.5 16];
    elseif strcmp(syl,'c')
        f_cutoff=[250 3000];
        t_assay=80;
        spect_params=[0 8];
    elseif strcmp(syl,'d')
        f_cutoff=[3000 5000];
        t_assay=.004;
        spect_params=[0 8];
    elseif strcmp(syl,'n')
        f_cutoff=[2000 5000];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'i')
        f_cutoff=[0 2000];
        t_assay=50;
        spect_params=[0 8];
    end
elseif strcmp(bname,'pu54b57')   % Jon's bird
    if strcmp(syl,'a')
        f_cutoff=[1000 3000];
        t_assay=50;
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[1500 2500];
        t_assay=.018;
        spect_params=[.5 16];
    elseif strcmp(syl,'c')
        f_cutoff=[1500 2500];
        t_assay=.018;
        spect_params=[.5 16];
    elseif strcmp(syl,'f')
        f_cutoff=[1500 2500];
        t_assay=.018;
        spect_params=[.5 16];
    elseif strcmp(syl,'d')
        f_cutoff=[3000 4500];
        t_assay=.024;
        spect_params=[0.5 16];
    elseif strcmp(syl,'j')
        f_cutoff=[3100 4500];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'n')
        f_cutoff=[3100 4500];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'i')
        f_cutoff=[3000 4500];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'e')
        f_cutoff=[1000 1650];
        t_assay=80;
        spect_params=[0.5 16];
    end
elseif strcmp(bname,'pu53b58')   % Jon's bird
    if strcmp(syl,'a')
        f_cutoff=[1000 2400];
        t_assay=50;
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[1500 2500];
        t_assay=.018;
        spect_params=[.5 16];
    elseif strcmp(syl,'j')
        f_cutoff=[3000 4500];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'n')
        f_cutoff=[3000 4500];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'k')
        f_cutoff=[800 1200];
        t_assay=50;
        spect_params=[0.5 16];
    end
elseif strcmp(bname,'pu28b75')   % Jon's bird
    if strcmp(syl,'a')
        f_cutoff=[500 1200];
        t_assay=65;
        spect_params=[0.5 16];
    elseif strcmp(syl,'b')
        f_cutoff=[3000 4200];
        t_assay=40;
        spect_params=[0.5 16];
    elseif strcmp(syl,'c')
        f_cutoff=[3000 4500];
        t_assay=.010;
        spect_params=[0.5 16];
    elseif strcmp(syl,'d')
        f_cutoff=[800 1200];
        t_assay=50;
        spect_params=[0.5 16];
    elseif strcmp(syl,'e')
        f_cutoff=[3000 5000];
        t_assay=.018;
        spect_params=[0.5 16];
    elseif strcmp(syl,'f')
        f_cutoff=[3000 5000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'x')
        f_cutoff=[500 1500];
        t_assay=40;
        spect_params=[0.5 16];
    end
elseif strcmp(bname,'ZFpi90w37')   % NCM impland
    if strcmp(syl,'a')
        f_cutoff=[2250 2850];
        t_assay=.024;
        spect_params=[0 16];
    elseif strcmp(syl,'b')
        f_cutoff=[2750 3750];
        t_assay=.04;
        spect_params=[0 16];
    elseif strcmp(syl,'c')
        f_cutoff=[4000 8000];
        t_assay=50;
        spect_params=[.5 16];
    elseif strcmp(syl,'d')
        f_cutoff=[2000 3500];
        t_assay=.008;
        spect_params=[0 16];
    elseif strcmp(syl,'e')
        f_cutoff=[4000 6000];
        t_assay=.008;
        spect_params=[0 16];
    elseif strcmp(syl,'f')
        f_cutoff=[4000 8000];
        t_assay=50;
        spect_params=[.5 16];
    elseif strcmp(syl,'g')
        f_cutoff=[2000 4000];
        t_assay=80;
        spect_params=[0 16];
    elseif strcmp(syl,'h')
        f_cutoff=[2000 2850];
        t_assay=.04;
        spect_params=[0 16];
    end
elseif strcmp(bname,'g44o23')
    if strcmp(syl,'a')
        if 0
            % original
            f_cutoff=[100 3000];
            t_assay=25;
            spect_params=[0 8];
        elseif 1
            % newer - trying to get a more consitent meaurement of
            % pitch at start of downsweep.  changed to this 10/24/2006
            f_cutoff=[1250 3000];
            t_assay=.008;
            spect_params=[0.5 16];
        end
    elseif strcmp(syl,'b')
        if 0
            % original
            f_cutoff=[3000 5000];
            t_assay=.012;
            spect_params=[0 8];

            % problem with above is existence of 2 features
            % at similar powers and about 1000 hz apart - 4k and 5k

        elseif 1
            % newest - go for high harmonic
            f_cutoff=[7000 9500];
            t_assay=.012;
            spect_params=[0 8];
        else
            % new - goes for initial, short, higher pitched element
            f_cutoff=[3500 5250];
            t_assay=.004;
            spect_params=[0 8];
        end
    elseif strcmp(syl,'c')
        f_cutoff=[3000 5000];
        t_assay=.016;
        spect_params=[.5 16];
    elseif strcmp(syl,'d')
        f_cutoff=[4000 8000];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'e')
        f_cutoff=[500 1750];
        t_assay=50;
        spect_params=[0 8];
    elseif strcmp(syl,'i')
        f_cutoff=[2000 5000];
        t_assay=50;
        spect_params=[0 8];
    elseif strcmp(syl,'j')
        f_cutoff=[750 1500];
        t_assay=50;
        spect_params=[0.5 16];
    end
elseif strcmp(bname,'pu26y2')
    if strcmp(syl,'b')
        f_cutoff=[1500 3500];
        t_assay=.016;
        spect_params=[0.5 16];
    elseif strcmp(syl,'c')
        f_cutoff=[2000 6000];
        t_assay=.016;
        spect_params=[0.5 16];
% the below might work for quantifying second half of syl
% f_cutoff=[500 2450];
%         t_assay=.07;
% spect_params=[0 8];
    elseif strcmp(syl,'d')
%        f_cutoff=[2000 4000];
% changed low limit 12/13/06 - a few low errors.
        f_cutoff=[2500 4000];
        t_assay=.02;
        spect_params=[0 8];
    elseif strcmp(syl,'e')
        f_cutoff=[500 9000];
        t_assay=.012;
        spect_params=[0 8];        
    elseif strcmp(syl,'f')
        % syl f is not exactly like 'b' - first part is noisier
        f_cutoff=[1500 3500];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'g')
        % syl g MIGHT be different from 'c' - not sure 
        f_cutoff=[2000 5000];
        t_assay=.024;
        spect_params=[0.5 16];
    elseif strcmp(syl,'j')
        f_cutoff=[500 2000];
        t_assay=.01;
        spect_params=[.5 16];
        % below doesnt really measure anything - a dummy var set
    elseif strcmp(syl,'a')
        f_cutoff=[2000 2500];
        t_assay=.004;
        spect_params=[0 8];
    end
elseif strcmp(bname,'r12r11')
    % a - very low pitch, not great resolution since few freq points
    % (although doubling resolution doesnt seem to change mean).
    %
    % b - range probably needs to be a bit higher to capture upshifted
    % song.  as is, would benefit a bit from only using within-range peaks
    % (not edge peaks).
    %
    % m - a LOT of edge peaks at top and bottom. definitely do
    % within-range peak analysis.
    %
    % n - edge peaks not a worry.  but this stack is sloped some of the
    % time, and slope seems to change with shifts.  worried that pitch-quant
    % time (which is late in syl) is bouncing around to different parts of
    % stack.
    %
    % d - extremely noisy and hard to quantify.  as is, trying to quant a very
    % low pitched element (few freq points).  does have a lot of edge
    % peaks; see if this helps.
    %
    % e - quant seems okrywrtyrtywr
    if strcmp(syl,'a')
        f_cutoff=[900 1300];
        t_assay=.016;
        spect_params=[0.5 16];

        % better resolution, but for 1_15 and 1_16 at least doesnt affect
        % mean very much
% f_cutoff=[900 1300];
%         t_assay=.016;
%         spect_params=[0.5 32];
    elseif strcmp(syl,'b')
        f_cutoff=[3000 4150];
        t_assay=.03;
        spect_params=[0.5 16];
    elseif strcmp(syl,'m')
        % for first stack in syllable 'c'
        f_cutoff=[2200 2760];
        t_assay=.07;
        spect_params=[0.5 16];
    elseif strcmp(syl,'n')
        % for second stack  in syllable 'c'
                f_cutoff=[5700 7000];
        t_assay=.12;
        spect_params=[0.5 16];
    elseif strcmp(syl,'d')
%         f_cutoff=[500 1750];
%         t_assay=.016;
%         spect_params=[0.5 16];

%        f_cutoff=[1000 1600];
        f_cutoff=[1100 1500];
        t_assay=.024;
        spect_params=[0.5 16];
    elseif strcmp(syl,'e')
        f_cutoff=[3000 4200];
        t_assay=.024;
        spect_params=[0.5 16];
    elseif strcmp(syl,'w')      % calls
        f_cutoff=[1000 8000];
        t_assay=.004;
        spect_params=[0 8];
    end
elseif strcmp(bname,'pk7r88')
    % a has some problems - pitch quant tries to hit high element in second
    % 8ms bin.  if first part of syl is truncated, this gets missed.
    % also, the high element sometimes bifurcates - a bit of a bimodal
    % dist.  a couple of scattered incidences of edge peaks, but these seem
    % to be in truncated syls
    %
    % b - looks ok, no edge effects.
    %
    % c - no trouble with edges. a stack with a small slope (could be an
    % issue).  looks good though
    %
    % d - clean harmonics but VERY low pitched.  a LOT of upper-edge peaks,
    % run with new method for sure.
    %
    % e - a fairly noisy syl but quant is ok.  a few outliers but not a big
    % problem.
    %
    % f - problems late in experiment - quieter/shortes syl, missing the
    % first part of a very sloped syl leads to a subset having too-low
    % pitch.  just get rid of quants between 1000-1600
    %
    % g - big problems with high edge peak - reduce upper limit to 1100 and
    % rerun% REDID PARAMS - now good
    if strcmp(syl,'a') 
%        f_cutoff=[2000 5400];
        f_cutoff=[3500 5400];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'b')
        f_cutoff=[1500 2500];
        t_assay=.03;
        spect_params=[0.5 16];
    elseif strcmp(syl,'c')
        f_cutoff=[5750 7400];
        t_assay=.016;
        spect_params=[0.5 16];
    elseif strcmp(syl,'d')  % tough one
        f_cutoff=[500 1005];
        t_assay=.016;
        spect_params=[0.5 16];
    elseif strcmp(syl,'e')
        f_cutoff=[1500 3000];
        t_assay=.016;
        spect_params=[0.5 16];
    elseif strcmp(syl,'f')
        f_cutoff=[500 3250];
        t_assay=.012;
        spect_params=[0 8];
    elseif strcmp(syl,'g')  
%         f_cutoff=[200 1100];
%         t_assay=60;
%         spect_params=[0.5 16];
        f_cutoff=[3000 5000];
        t_assay=.012;
        spect_params=[0 8];
    end
end



if ~exist('f_cutoff')    % if undefined
    f_cutoff='undefined';
    t_assay='undefined';
    spect_params='undefined';
end

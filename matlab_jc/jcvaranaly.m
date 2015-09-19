% Figure 1 sensitivity analysis  %%%%%%%%%%%%%%%%%%%%%
        % A. Note onset time (shift right or left 5ms - note that the 0 value is
        %  1.5*(minimum sd) of sd curve.
            % Short notes are fairly sensitive to this. Long notes are very very robust.
        % B. x-axis as % of note vs. distance in ms from front of note
            % Both types of notes are robust.  Short notes prefer ms from front
            % whereas long notes prefer % of note.  Stick with % of note for
            % further analysis.
        % C. Definition of shifted - prct of max shift vs. time from center of max
        % shift dataset
        %   High robustness to this. Stick with 80% of max shift.
        % D. divide by std curve vs. don't - similar results - actually higher
        % correlations for both without doing it!!!
        [bigtarget,bigshift]=targetingplot(Alldata,1,0.8,2,0);
        figure;xlim([0 1]);ylim([0 1]);plot(bigtarget,bigshift,'*')
        hold on;plot(bigtarget(Alldata(1).ind_longnotes),bigshift(Alldata(1).ind_longnotes),'*','Color','red')
        corrlong=corrcoef(bigtarget(Alldata(1).ind_longnotes),bigshift(Alldata(1).ind_longnotes))
        corrshort=corrcoef(bigtarget(Alldata(1).ind_shortnotes),bigshift(Alldata(1).ind_shortnotes))
        corrs=corrcoef(bigtarget,bigshift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 2 %%%%%%%%%%%%%%%%%%%%%
        figure;targetingplotFigure2(Alldata2(1).ind_longnotes,Alldata2,1,0.8,2,0);
        % old method - [avgpcurvetrace]=targetingplotfig2(Alldata,1,0.8,2,0);

        % is not more pronounced earlier on in the
        % shifting protocol (e.g. when it hits 1sd above)

        % One might have expected a strong inverse relationship between std of
        % targeting and precision of pitch bump, but only weak, non-significant
        % relationships were observed and these relationships were inconsistent
        % across different ways of measuring precision of pitch bump.  
        % Another possible explanation for the variation in pitch bump precision is variation
        % in the time course of baseline variation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Figure 3 and sensitivity analysis %%%%%%%%%%%%%%%%%%%%

        %%%%% Print curves at one contingency above/below (30=20 and 80)
        Figure3analy(3,Alldata,30);
        xlim([0 25]);ylim([0 1.05])

        %%% Test sensitivity to different contingencies
        [BSlopesLong,BSlopesShort]=Figure3Banaly(Alldata,[0 5 10 15 30 25 30 35 40]);
        figure;plot(BSlopesShort(:,1)/BSlopesLong(:,1),'k')
        hold on;plot(BSlopesShort(:,2)/BSlopesLong(:,2),'r')
        hold on;plot(BSlopesShort(:,3)/BSlopesLong(:,3),'b')
        % or
        figure;plot(BSlopesShort(:,1),'k')
        hold on;plot(BSlopesLong(:,1),'k')
        hold on;plot(BSlopesShort(:,2),'r')
        hold on;plot(BSlopesLong(:,2),'r')
        hold on;plot(BSlopesShort(:,3),'b')
        hold on;plot(BSlopesLong(:,3),'b')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

Alldata(1).exp=BK15BK14_exp1;
Alldata(2).exp=BK15BK14_exp2;
Alldata(3).exp=BK15BK14_exp3;
Alldata(4).exp=BK15BK14_exp4;
Alldata(5).exp=BK20BK45_A;
Alldata(6).exp=BK20BK45_B;
Alldata(7).exp=BK50W18_endup;
Alldata(8).exp=BK50W18_frontdown;
Alldata(9).exp=BK50W18_frontsup;
Alldata(10).exp=PK16R47;
Alldata(11).exp=PK20R49;
Alldata(12).exp=PK32BK28_exp2;
Alldata(13).exp=PK37BK19_middledown;
Alldata(14).exp=PK37BK19_middleup;
Alldata(15).exp=PK39BK5;
Alldata(16).exp=PU34;
Alldata(17).exp=OR92OR82;

% Get all baseline acsf and inactivation data into files
figure;targetingplotFigure2(Alldata2(1).ind_longnotes,Alldata2,1,0.8,2,0);
Alldata=baselinecat(Alldata);

for i=1:length(Alldata)
    figure;plot(std(Alldata(i).baselineAC'))
end

startnote=[668 800 503 462 309 228 221 442 392 435 239 246 254 273 257 218];
endnote=[1305 1393 1088 1008 971 440 455 866 815 867 852 419 480 793 443 852];
for i=1:length(Alldata)
    Alldata(i).startnote=startnote(17-i);
    Alldata(i).endnote=endnote(17-i);
end
for i=1:length(Alldata)
    ons=Alldata(i).startnote;
    offs=Alldata(i).endnote;
    exp=Alldata(i).exp;
    count=0;
    x=[];
    y=[];
    for j=1:length(exp)
        if exp(1).acsf(j)==1
            count=count+1;
            x(count)=(exp(1).tfromwnonend(j)+exp(1).tfromwnonbegin(j))/2;
            y(count)=mean(mean(exp(j).selectedpitchcurves(ons+80:offs-80,:)'));
        end
    end
    figure;plot(x,y)
end

notshift=0;
shiftedmethod=1;
prct=0.1;
for i=1:9
    [bigtarget,bigshift]=targetingplot(Alldata,shiftedmethod,prct,0,notshift);
    c=corrcoef(bigtarget,bigshift);
    corrmat(i)=c(2);
    prct=prct+-0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;targetingplotFigure2(Alldata2(1).ind_longnotes,Alldata2,1,0.8,2,0);
%get toffs
n=16;
toffs=[];
for i=7:20
    toffs=[toffs;Alldata2(n).exp(i).toffset];
end
figure;plot(toffs)
%
toffs1=[];
for j=1:length(toffs)
    if toffs(j)>0
        toffs1(j)=toffs(j);
    else
        toffs1(j)=mean(toffs);
    end
end

toffsets(n).data=toffs;

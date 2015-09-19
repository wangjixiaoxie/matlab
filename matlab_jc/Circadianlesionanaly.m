% 1. For all data pre/post - p=0.053 for destabilization post-lesion,
% where destabilization is defined as abs(Slope_post)-abs(Slope_pre)>0
% There is no trend for the slope to become positive or negative - both
% directions of destabilization are equally likely.  

% 2A. Breaking the data up into post_early and post_late does not reveal a
% dichotomy - if anything, the effect is slightly more pronounced
% post_early.

% 2B. There is actually a trend toward stabilization over time from pre_early to
% pre_late (p=0.15) and a slight trend toward stabilization over time from post_early to
% post_late (p=0.33).  So it's quite surprising to see such a
% destabilizing effect right after the lesion.  But it doesn't increase the
% significance of things to look pre_late vs. post_early, presumably
% because there is much less data and more variation.

% 4. Looking at the actual data to see if destabilization is linear...


% bn(1).bird=o50bk72;
% bn(2).bird=bk72w64;
% bn(3).bird=bk13w63;
% bn(4).bird=bk35w47;
% bn(5).bird=pk20bk78;
% bn(6).bird=pk28k66;
% bn(7).bird=pk29bk36;
% bn(8).bird=pk42bk73;
% bn(9).bird=pk30bk79;
% bn(10).bird=o7bk77;
% bn(11).bird=pk80r6;

% Simple linear regression - How much does each hour of time tell you about FF in Hz?
for i=1:length(bn)
    % convert to hours - done
    % calculate sums
    Sxpre=sum(tvalspres(i).data);Sypre=sum(Fvalspres(i).data);Sxxpre=sum(tvalspres(i).data.*tvalspres(i).data);
    Sxypre=sum(tvalspres(i).data.*Fvalspres(i).data);Syypre=sum(Fvalspres(i).data.*Fvalspres(i).data);
    Sxpost=sum(tvalsposts(i).data);Sypost=sum(Fvalsposts(i).data);Sxxpost=sum(tvalsposts(i).data.*tvalsposts(i).data);
    Sxypost=sum(tvalsposts(i).data.*Fvalsposts(i).data);Syypost=sum(Fvalsposts(i).data.*Fvalsposts(i).data);
    npre=length(tvalspres(i).data);npost=length(tvalsposts(i).data);
    % calculate slope estimates
    Beta_pre=(npre*Sxypre-Sxpre*Sypre)/(npre*Sxxpre-Sxpre*Sxpre);
    Beta_post=(npost*Sxypost-Sxpost*Sypost)/(npost*Sxxpost-Sxpost*Sxpost);
    % calculate standard errors
    seE_pre2=(1/(npre*(npre-2)))*(npre*Syypre-Sypre*Sypre-Beta_pre*Beta_pre*(npre*Sxxpre-Sxpre*Sxpre));
    seE_post2=(1/(npost*(npost-2)))*(npost*Syypost-Sypost*Sypost-Beta_post*Beta_post*(npost*Sxxpost-Sxpost*Sxpost));
    seB_pre(i)=sqrt((npre*seE_pre2)/(npre*Sxxpre-Sxpre*Sxpre));
    seB_post(i)=sqrt((npost*seE_post2)/(npost*Sxxpost-Sxpost*Sxpost));
    Beta_pres(i)=Beta_pre;
    Beta_posts(i)=Beta_post;
    if abs(Beta_post)>abs(Beta_pre)
        evidenceforSLR(i)=(Beta_post-Beta_pre)/mean(seB_pre(i)+seB_post(i)); % t-scores
    else
        evidenceagainstSLR(i)=(Beta_pre-Beta_post)/mean(seB_pre(i)+seB_post(i)); % t-scores
    end
end

%%% Linear covariance approach - what I did originally - Gives the exact same answer
    % This demonstrates that amount of variation in FF doesn't matter
for i=1:length(bn)
        Cmatpre=cov(tvalspres(i).data,Fvalspres(i).data);
        % covariance of F and t, divided by covariance of t and t (i.e. variance of t) - 
            % how much F varies with t normalized by how much t varies (i.e. rise/run)
        slopepre(i)=Cmatpre(2)/Cmatpre(1); 
        Cmatpost=cov(tvalsposts(i).data,Fvalsposts(i).data);
        slopepost(i)=Cmatpost(2)/Cmatpost(1);
end
    % This approach is robust to differences in variation pre vs. post - mean(mm)=2, which is the true trend
                %         for i=1:500
                %             F1=rand(1,1000)*1000+[2:2:2000];
                %             F2=[1:1:1000];
                %             F1=F1-mean(F1);
                %             F2=F2-mean(F2);
                %             mm(i)=mean((F1.*F2)./(F2.*F2));
                %         end
    % This approach allows you to make 95% confidence intervals on pre and
    % then determine if post falls within these intervals:
            timeall=[timepre;timepost];
            freqall=[FFpre;FFpost];
%             proppost=1-proppre;
                for i=1:10000
                    % sample from pre
                    randpre=[ceil(rand(1,500)*(length(timepre)))];
                    % sample from post
                    randpost=[floor(length(timepre)+rand(1,500)*length(timepost))];
                    randsamp=[randpre randpost];
                    FFtest=freqall(randsamp)-mean(freqall(randsamp));
                    timetest=timeall(randsamp)-mean(timeall(randsamp));
                    slopetest(i)=sum(timetest.*FFtest)./sum(timetest.*timetest);
                end
                evidencefor(iii)=0;
                evidenceagainst(iii)=0;
                if abs(slope_post)>abs(slope_pre)
                    if slope_post>0
                        evidencefor(iii)=1-sum(slopetest>slope_post)/10000-1/10000;
                    else
                        evidencefor(iii)=1-sum(slopetest<slope_post)/10000-1/10000;
                    end
                else
                    if slope_pre>0
                        evidenceagainst(iii)=1-sum(slopetest>slope_pre)/10000-1/10000;
                    else
                        evidenceagainst(iii)=1-sum(slopetest<slope_pre)/10000-1/10000;
                    end
                end
                sloper_pre(iii)=slope_pre;
                sloper_post(iii)=slope_post;
end
                

bn(1).bird=o50bk72;
bn(2).bird=bk72w64;
bn(3).bird=bk13w63;
bn(4).bird=bk35w47;
bn(5).bird=pk20bk78;
bn(6).bird=pk28k66;
bn(7).bird=pk29bk36;
bn(8).bird=pk42bk73;
bn(9).bird=pk30bk79;
bn(10).bird=o7bk77;
bn(11).bird=pk80r6;
for iii=1:length(bn)
thisbird=bn(iii).bird;
figure;hold on;
subplot(211);hold on;
dpl=[];
dpl=zeros(1,length(thisbird));
for i=1:length(thisbird)
    if ~isempty(thisbird(i).dayspostlesion)
    dpl(i)=thisbird(i).dayspostlesion;
    end
end
preend=max(find(dpl<0));
element=1;
Fvalspre=[];
tvalspre=[];
while element<preend;
    element=element+1;
    if thisbird(element).keep 
        element
        plot(thisbird(element).timevals,thisbird(element).FFvals,'*')
        Fvalspre=[Fvalspre;thisbird(element).FFvals];
        tvalspre=[tvalspre;thisbird(element).timevals];
    end
end
element=element+1;
xlim([0 900])
% lowval=min(thisbird(count).FFvals)-100;
% highval=max(thisbird(count).FFvals)+100;
% ylim([lowval highval])
subplot(212);hold on;
Fvalspost=[];
tvalspost=[];
% for i=element:length(thisbird)
%     if thisbird(i).keep
%         plot(thisbird(i).timevals,thisbird(i).FFvals,'*')
%         Fvalspost=[Fvalspost;thisbird(i).FFvals];
%         tvalspost=[tvalspost;thisbird(i).timevals];
% 
%     end
% end
i=element;
while i<length(thisbird)
     if thisbird(i).keep %%&& thisbird(i).dayspostlesion<14
         plot(thisbird(i).timevals,thisbird(i).FFvals,'*')
        Fvalspost=[Fvalspost;thisbird(i).FFvals];
        tvalspost=[tvalspost;thisbird(i).timevals];
     end
    i=i+1;
end
        if iii==8
            Fvalspost=Fvalspost(700:end); % this is screwy;
            tvalspost=tvalspost(700:end);
        end
xlim([0 900])
% ylim([lowval highval])

emptiness(iii)=isempty(Fvalspost);


%%% Linear covariance approach - How much does each hour of time tell you about FF in Hz?
                FFpre=Fvalspre-mean(Fvalspre);
        timepre=tvalspre-mean(tvalspre); 
        timepre=timepre/60; % convert to hours
        slope_pre=sum(timepre.*FFpre)./sum(timepre.*timepre)  % tells you 3.34Hz pre
        FFpost=Fvalspost-mean(Fvalspost);
        timepost=tvalspost-mean(tvalspost); 
        timepost=timepost/60; % convert to hours
        slope_post=sum(timepost.*FFpost)./sum(timepost.*timepost) % tells you 6.88Hz post

    % This approach is robust to differences in variation pre vs. post - mean(mm)=2, which is the true trend
                %         for i=1:500
                %             F1=rand(1,1000)*1000+[2:2:2000];
                %             F2=[1:1:1000];
                %             F1=F1-mean(F1);
                %             F2=F2-mean(F2);
                %             mm(i)=mean((F1.*F2)./(F2.*F2));
                %         end
    % This approach allows you to make 95% confidence intervals on pre and
    % then determine if post falls within these intervals:
            timeall=[timepre;timepost];
            freqall=[FFpre;FFpost];
%             proppost=1-proppre;
                for i=1:10000
                    % sample from pre
                    randpre=[ceil(rand(1,500)*(length(timepre)))];
                    % sample from post
                    randpost=[floor(length(timepre)+rand(1,500)*length(timepost))];
                    randsamp=[randpre randpost];
                    FFtest=freqall(randsamp)-mean(freqall(randsamp));
                    timetest=timeall(randsamp)-mean(timeall(randsamp));
                    slopetest(i)=sum(timetest.*FFtest)./sum(timetest.*timetest);
                end
                evidencefor(iii)=0;
                evidenceagainst(iii)=0;
                if abs(slope_post)>abs(slope_pre)
                    if slope_post>0
                        evidencefor(iii)=1-sum(slopetest>slope_post)/10000-1/10000;
                    else
                        evidencefor(iii)=1-sum(slopetest<slope_post)/10000-1/10000;
                    end
                else
                    if slope_pre>0
                        evidenceagainst(iii)=1-sum(slopetest>slope_pre)/10000-1/10000;
                    else
                        evidenceagainst(iii)=1-sum(slopetest<slope_pre)/10000-1/10000;
                    end
                end
                sloper_pre(iii)=slope_pre;
                sloper_post(iii)=slope_post;
end
                
%             prctile((slopetest),97.5) % compare this to postval
%             prctile((slopetest),2.5) % compare this to postval
%     slope_pre
%     slope_post

        
            
            
            
            
            
%             
%             
% % Hypothesis 1: There is a circadian trend after lesions.
%     % Fit with first order polynomial (i.e. circadian trend) plus gaussian noise.
%     % p_post(1) vs p_pre(1)
%         [p_post,s_post,mu_post]=polyfit(tvalspost,Fvalspost,1)
%         [p_pre,s_pre,mu_pre]=polyfit(tvalspre,Fvalspre,1)
%         [y_post,delta_post]=polyval(p_post,tvalspost,s_post,mu_post);
%         [y_pre,delta_pre]=polyval(p_pre,tvalspre,s_pre,mu_pre);
%        % Stdev should be much smaller post after factoring in the trend
%         std(abs(Fvalspost-y_post))/std(abs(Fvalspre-y_pre))
%         std(Fvalspost)/std(Fvalspre)
%     % Presumably p shouldn't depend on CV because gaussian noise is included in the model?  Is this accurate?
%     
% % Hypothesis 2: 
% 
% 

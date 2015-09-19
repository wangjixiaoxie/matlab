function [outmat,ss]=learn_plot(sqout,ss)



    
    %first plot example cbin

    subplot(3,4,1:3)
%     ax=gca();
%     exsong.ax=ax;
%     exsong.path='/oriole7/dir2/pu19bk81/screen'
%     exsong.fn='pu19bk81_280311_1254.-29188.cbin';
%     exsong.bnds=[3.75 5.65]
%     plotcbin(exsong);
% 
    %this plots example transition probabilities for 6 days for this example
    subplot(3,4,5:7);
    ssind=2;
    
    cmd=['cd ' ss(ssind).pth]
    eval(cmd);
    cmd2='load sumdata.mat'
    eval(cmd2);
    inds=[4:11]
    trgnt=avls.SEQTRGNT{1};
    
    st.prob(:,1)=outnotect(inds,trgnt)./sum(outnotect(inds,:),2);
    st.prob(:,2)=1-st.prob(:,1);
    st.col{1}='k'
    st.col{2}='r'
    
    st.days=datenum(dayout(inds))-datenum(dayout(inds(1)))+1;
    st.err=[];
    plotseqtm(st);
    
    %load the repeat data
   
    
    subplot(3,4,8);
    clear st
    inds=[11:14]
    trgnt=avls.SEQTRGNT{1};
    
    st.prob(:,1)=outnotect(inds,trgnt)./sum(outnotect(inds,:),2);
    st.prob(:,2)=1-st.prob(:,1);
    st.col{1}='k'
    st.col{2}='r'
    
    st.days=datenum(dayout(inds))-datenum(dayout(inds(1)))+1;
    st.err=[];
    plotseqtm(st);
    
    
     cd ~/matlab/seqpap/learn_anal
    load transitionTW.mat
    rpout=transitionTW
    subplot(3,4,9:10)
  
    plotscatterlearn(sqout,rpout )
    subplot(3,4,11:12)
    plotscatterrec(sqout,rpout )
    
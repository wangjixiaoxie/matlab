% b10o7
% surgery on 6.25.09
% 7.23 - probes in
% 7.25 (morning) - started singing
% EXPERIMENT 1
    % 7.27 at 11am - 2mM AP5 at 1.0uL/min
    % 7.27 at 12:50pm - ampon - hit below 2120Hz 
    % 7.27 at 1pm - hit below 2150Hz
    % 7.27 at 5:02pm - door open - acquisition off - ACSF at 1.5uL/min
    % 7.27 at 5:32pm - door closed - acquisition on - ACSF at 0.6uL/min
    
    % Dose of AP5 only lowered CV by 15%, although it was fairly low to
    % begin with.  There was definite learning after AP5 went off, but
    % there also was some learning with AP5 present.  Thus I will try a
    % higher dose of AP5 for the next run.
% EXPERIMENT 2
    % 7.28 at 3:00pm - 4mM AP5 at 1.0uL/min
    % 7.28 at 4:52pm - ampon - hit below 2230Hz
    % 7.28 at 8:29pm - door open - acquisition off - ACSF at 1.5uL/min
    % 7.28 at 8:59pm - door closed - acquisition on - ACSFsmb://chick/Djcharles/b10o7/729_ACSFampoff/b10o7_290709_2033.739.rec at 0.6uL/min
    
% EXPERIMENT 3
    % 7.29 at 1:15pm - 4mM AP5 at 1.0uL/min
    % 7.29 at 3pm - new template
    % 7.29 at 3:16pm - ampon - hit below 2250Hz
    % 7.29 at 5:38pm - hit below 2270Hz
    % 7.29 at 6:48pm - door open - acquisition off - ACSF at 1.5uL/min
    % 7.29 at 7:18pm - door closed - acquisition on - ACSF at 0.6uL/min
    % 7.30 morning - clogged
    % 7.30 ~12pm - unplugged and handled to put epoxy on ends of probe
        % note that probes remain in brain to ctl for damage
figure;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv tvals729apvWN tvals729ampoff],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:)) median(pitch729apv(150:170,:)) median(pitch729apvWN(150:170,:)) median(pitch729ampoff(150:170,:))],'*','Color','b')         
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv tvals729apvWN],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:)) median(pitch729apv(150:170,:)) median(pitch729apvWN(150:170,:))],'*','Color','g')     
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:)) median(pitch729apv(150:170,:))],'*','Color','r')   
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:))],'*','Color','b')   
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:))],'*','Color','b')

valsX=evtaf_freq('batch',[1800,2500],'a',128,'obs0',1,1);


%%%%%%%%
%%%%%%% Figure
figure;subplot(311);hold on;
plot(median(pitch727A'))
hold on;plot(median(pitch728off'),'k')
hold on;plot(median(pitch728A'),'r')
hold on;plot(median(pitch729ampoff'),'g')
subplot(312);hold on;
plot(median(pitch727A'))
hold on;plot(median(pitch728off'),'k')
hold on;plot(median(pitch728A'),'r')
hold on;plot(median(pitch729ampoff'),'g')
plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv tvals729apvWN tvals729ampoff],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:)) median(pitch729apv(150:170,:)) median(pitch729apvWN(150:170,:)) median(pitch729ampoff(150:170,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv tvals729apvWN],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:)) median(pitch729apv(150:170,:)) median(pitch729apvWN(150:170,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:)) median(pitch729apv(150:170,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:)) median(pitch728off(150:170,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:)) median(pitch728apvAmpon(150:170,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:)) median(pitch728apv(150:170,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:)) median(pitch727ampoff(150:170,:)) median(pitch728A(150:170,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:)) median(pitch727hits(150:170,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:)) median(pitch727APV(150:170,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A],[median(pitch725(150:170,:)) median(pitch726A(150:170,:)) median(pitch726B(150:170,:)) median(pitch727A(150:170,:))],'*','Color','b')
subplot(313);hold on;
plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv tvals729apvWN tvals729ampoff],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:)) median(pitch728apv(1600:1620,:)) median(pitch728apvAmpon(1600:1620,:)) median(pitch728off(1600:1620,:)) median(pitch729apv(1600:1620,:)) median(pitch729apvWN(1600:1620,:)) median(pitch729ampoff(1600:1620,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv tvals729apvWN],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:)) median(pitch728apv(1600:1620,:)) median(pitch728apvAmpon(1600:1620,:)) median(pitch728off(1600:1620,:)) median(pitch729apv(1600:1620,:)) median(pitch729apvWN(1600:1620,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off tvals729apv],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:)) median(pitch728apv(1600:1620,:)) median(pitch728apvAmpon(1600:1620,:)) median(pitch728off(1600:1620,:)) median(pitch729apv(1600:1620,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon tvals728off],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:)) median(pitch728apv(1600:1620,:)) median(pitch728apvAmpon(1600:1620,:)) median(pitch728off(1600:1620,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv tvals728apvAmpon],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:)) median(pitch728apv(1600:1620,:)) median(pitch728apvAmpon(1600:1620,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A tvals728apv],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:)) median(pitch728apv(1600:1620,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits tvals727ampoff tvals728A],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:)) median(pitch727ampoff(1600:1620,:)) median(pitch728A(1600:1620,:))],'*','Color','b')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV tvals727hits],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:)) median(pitch727hits(1600:1620,:))],'*','Color','g')
hold on;plot([tvals725 tvals726A tvals726B tvals727A tvals727APV],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:)) median(pitch727APV(1600:1620,:))],'*','Color','r')
hold on;plot([tvals725 tvals726A tvals726B tvals727A],[median(pitch725(1600:1620,:)) median(pitch726A(1600:1620,:)) median(pitch726B(1600:1620,:)) median(pitch727A(1600:1620,:))],'*','Color','b')
%%%%%%%%%%%
%%%%%%%%%%%
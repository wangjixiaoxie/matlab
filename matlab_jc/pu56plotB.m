% Experiment 5 - control - pu56
    pitchPreExp5=[pitch510,pitch511A];
    pitchWNExp5=[pitch511WN];
    pitchPostExp5=[pitch514wnoff];
    timevalsPreExp5=[timevals510,timevals511A];
    timevalsWNExp5=[timevals511WN];
    timevalsPostExp5=[timevals514wnoff];
    figure;plot([timevalsPreExp5 timevalsWNExp5 timevalsPostExp5],[median(pitchPreExp5(160:200,:)),median(pitchWNExp5(160:200,:)),median(pitchPostExp5(160:200,:))],'*','Color','k')
    hold on;plot([timevalsPreExp5 timevalsWNExp5],[median(pitchPreExp5(160:200,:)),median(pitchWNExp5(160:200,:))],'*','Color','r')
    plot([timevalsPreExp5],[median(pitchPreExp5(160:200,:))],'*','Color','k')
    wc=40;
    for i=1:floor(size(pitchPreExp5,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPreExp5(160:200,start:last)))];
        plot([timevalsPreExp5(start) timevalsPreExp5(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchWNExp5,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchWNExp5(160:200,start:last)))];
        plot([timevalsWNExp5(start) timevalsWNExp5(last)],[val1 val1],'Color','k','LineWidth',8)
    end
    for i=1:floor(size(pitchPostExp5,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPostExp5(160:200,start:last)))];
        plot([timevalsPostExp5(start) timevalsPostExp5(last)],[val1 val1],'Color','r','LineWidth',8)
    end

    % Contingencies
    hold on;plot([timevalsWNExp5(1) timevalsWNExp5(end)],[3700 3700],'Color','g','LineWidth',4)


% Experiment 4 - pu56 - 0.6uM - upshift
    pitchPreExp4=[pitch504A,pitch505A];
    pitchAPVExp4=[pitch505APV];
    pitchAPVwnExp4=[pitch505wn,pitch505wn2,pitch505wn3];
    pitchPostExp4=[pitch506wnoff,pitch510];
    timevalsPreExp4=[timevals504A,timevals505A];
    timevalsAPVExp4=[timevals505APV];
    timevalsAPVwnExp4=[timevals505wn,timevals505wn2,timevals505wn3];
    timevalsPostExp4=[timevals506wnoff,timevals510];
    figure;plot([timevalsPreExp4 timevalsAPVExp4 timevalsAPVwnExp4 timevalsPostExp4],[median(pitchPreExp4(160:200,:)),median(pitchAPVExp4(160:200,:)),median(pitchAPVwnExp4(160:200,:)),median(pitchPostExp4(160:200,:))],'*','Color','k')
    hold on;plot([timevalsPreExp4 timevalsAPVExp4 timevalsAPVwnExp4],[median(pitchPreExp4(160:200,:)),median(pitchAPVExp4(160:200,:)),median(pitchAPVwnExp4(160:200,:))],'*','Color','r')
    plot([timevalsPreExp4 timevalsAPVExp4],[median(pitchPreExp4(160:200,:)),median(pitchAPVExp4(160:200,:))],'*')
    plot([timevalsPreExp4],[median(pitchPreExp4(160:200,:))],'*','Color','k')

    wc=20;
    for i=1:floor(size(pitchPreExp4,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPreExp4(160:200,start:last)))];
        plot([timevalsPreExp4(start) timevalsPreExp4(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchAPVExp4,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchAPVExp4(160:200,start:last)))];
        plot([timevalsAPVExp4(start) timevalsAPVExp4(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchAPVwnExp4,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchAPVwnExp4(160:200,start:last)))];
        plot([timevalsAPVwnExp4(start) timevalsAPVwnExp4(last)],[val1 val1],'Color','k','LineWidth',8)
    end
    for i=1:floor(size(pitchPostExp4,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPostExp4(160:200,start:last)))];
        plot([timevalsPostExp4(start) timevalsPostExp4(last)],[val1 val1],'Color','r','LineWidth',8)
    end

    % Contingencies
    hold on;plot([timevals505wn(1) timevals505wn(end)],[3690 3690],'Color','g','LineWidth',8)
    hold on;plot([timevals505wn2(1) timevals505wn2(end)],[3670 3670],'Color','g','LineWidth',8)
    hold on;plot([timevals505wn3(1) timevals505wn3(end)],[3750 3750],'Color','g','LineWidth',8)




% Experiment 3 - pu56 - 1mM - downshift
    
    pitchPreExp3=[pitchPostExp2];
    timevalsPreExp3=[timevalsPostExp2];
    pitchPostExp3=[pitch501ac,pitch503A,pitch504A];
    timevalsPostExp3=[timevals501ac timevals503A timevals504A];
    figure;plot([timevalsPreExp3 timevals430apv timevals430apvwn timevalsPostExp3],[median(pitchPreExp3(160:200,:)) median(pitch430apv(160:200,:)) median(pitch430apvwn(160:200,:)) median(pitchPostExp3(160:200,:))],'*','Color','k')
    hold on;plot([timevalsPreExp3 timevals430apv timevals430apvwn],[median(pitchPreExp3(160:200,:)) median(pitch430apv(160:200,:)) median(pitch430apvwn(160:200,:))],'*','Color','r')
    plot([timevalsPreExp3 timevals430apv],[median(pitchPreExp3(160:200,:)) median(pitch430apv(160:200,:))],'*')
    plot([timevalsPreExp3],[median(pitchPreExp3(160:200,:))],'*','Color','k')
    
    wc=40;
    for i=1:floor(size(pitchPreExp3,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPreExp3(160:200,start:last)))];
        plot([timevalsPreExp3(start) timevalsPreExp3(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitch430apv,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitch430apv(160:200,start:last)))];
        plot([timevals430apv(start) timevals430apv(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitch430apvwn,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitch430apvwn(160:200,start:last)))];
        plot([timevals430apvwn(start) timevals430apvwn(last)],[val1 val1],'Color','k','LineWidth',8)
    end
    for i=1:floor(size(pitchPostExp3,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPostExp3(160:200,start:last)))];
        plot([timevalsPostExp3(start) timevalsPostExp3(last)],[val1 val1],'Color','r','LineWidth',8)
    end

    % Contingencies
    hold on;plot([timevals430apvwn(1) timevals430apvwn(37)],[3550 3550],'Color','g','LineWidth',4)
    hold on;plot([timevals430apvwn(38) timevals430apvwn(61)],[3650 3650],'Color','g','LineWidth',4)
    hold on;plot([timevals430apvwn(62) timevals430apvwn(end)],[3600 3600],'Color','g','LineWidth',4)

    
% Experiment 2 - pu56 - 1mM - upshift
    pitchPreExp2=[pitchPost,pitch1mMAPV,pitchAC426];
    pitchPostExp2=[pitch428,pitch428acwnoff,pitch430ac];
    timevalsPostExp2=[timevals428 timevals428acwnoff timevals430ac];
    timevalsPreExp2=[timevalsPost timevals1mMAPV timevalsAC426];
    
    figure;plot([timevalsPost timevals1mMAPV timevalsAC426 timevalsAPV427 timevalsAPVwn427 timevals428 timevals428acwnoff timevals430ac],[median(pitchPost(160:200,:)) median(pitch1mMAPV(160:200,:)) median(pitchAC426(160:200,:)) median(pitchAPV427(160:200,:)) median(pitchAPVwn427(160:200,:)) median(pitch428(160:200,:)) median(pitch428acwnoff(160:200,:)) median(pitch430ac(160:200,:))],'*','Color','k')
    hold on;plot([timevalsPost timevals1mMAPV timevalsAC426 timevalsAPV427 timevalsAPVwn427],[median(pitchPost(160:200,:)) median(pitch1mMAPV(160:200,:)) median(pitchAC426(160:200,:)) median(pitchAPV427(160:200,:)) median(pitchAPVwn427(160:200,:))],'*','Color','r')
    plot([timevalsPreExp2 timevalsAPV427],[median(pitchPost(160:200,:)) median(pitch1mMAPV(160:200,:)) median(pitchAC426(160:200,:)) median(pitchAPV427(160:200,:))],'*')
    plot([timevalsPreExp2],[median(pitchPreExp2(160:200,:))],'*','Color','k')
    
    wc=25;
    for i=1:floor(size(pitchPreExp2,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPreExp2(160:200,start:last)))];
        plot([timevalsPreExp2(start) timevalsPreExp2(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchAPV427,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchAPV427(160:200,start:last)))];
        plot([timevalsAPV427(start) timevalsAPV427(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchAPVwn427,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchAPVwn427(160:200,start:last)))];
        plot([timevalsAPVwn427(start) timevalsAPVwn427(last)],[val1 val1],'Color','k','LineWidth',8)
    end
    for i=1:floor(size(pitchPostExp2,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPostExp2(160:200,start:last)))];
        plot([timevalsPostExp2(start) timevalsPostExp2(last)],[val1 val1],'Color','r','LineWidth',8)
    end

    % Contingencies
    hold on;plot([timevalsAPVwn427(1) timevalsAPVwn427(102)],[3650 3650],'Color','g','LineWidth',4)
    hold on;plot([timevalsAPVwn427(103) timevalsAPVwn427(106)],[3600 3600],'Color','g','LineWidth',4)
    hold on;plot([timevalsAPVwn427(107) timevalsAPVwn427(end)],[3670 3670],'Color','g','LineWidth',4)

% Experiment 1- pu56 - 2mM - downshift
    figure;plot([timevalsProbein timevalsAPV timevalsAPVwnon timevalsPost],[median(pitchProbein(160:200,:)) median(pitchAPV(160:200,:)) median(pitchAPVwnon(160:200,:)) median(pitchPost(160:200,:))],'*','Color','k')
    hold on;plot([timevalsProbein timevalsAPV timevalsAPVwnon],[median(pitchProbein(160:200,:)) median(pitchAPV(160:200,:)) median(pitchAPVwnon(160:200,:))],'*','Color','r')
    hold on;plot([timevalsProbein timevalsAPV],[median(pitchProbein(160:200,:)) median(pitchAPV(160:200,:))],'*','Color','b')
    hold on;plot([timevalsProbein],[median(pitchProbein(160:200,:))],'*','Color','k')
% 
%     medpitchExp1=[median(pitchProbein(160:200,:)) median(pitchAPV(160:200,:)) median(pitchAPVwnon(160:200,:)) median(pitchPost(160:200,:))];
%     timevalsExp1=[timevalsProbein timevalsAPV timevalsAPVwnon timevalsPost];

    for i=1:floor(size(pitchProbein,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchProbein(160:200,start:last)))];
        plot([timevalsProbein(start) timevalsProbein(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchAPV,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchAPV(160:200,start:last)))];
        plot([timevalsAPV(start) timevalsAPV(last)],[val1 val1],'Color','r','LineWidth',8)
    end
    for i=1:floor(size(pitchAPVwnon,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchAPVwnon(160:200,start:last)))];
        plot([timevalsAPVwnon(start) timevalsAPVwnon(last)],[val1 val1],'Color','k','LineWidth',8)
    end
    for i=1:floor(size(pitchPost,2)/wc)
        start=i*wc-wc+1;
        last=start+wc-1;
        val1=[median(median(pitchPost(160:200,start:last)))];
        plot([timevalsPost(start) timevalsPost(last)],[val1 val1],'Color','r','LineWidth',8)
    end[
    % Contingencies
    hold on;plot([timevalsAPVwnon(1) timevalsAPVwnon(end)],[3550 3550],'Color','g','LineWidth',4)

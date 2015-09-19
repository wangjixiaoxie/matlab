% 1. BK15BK14 experiment 1 (April 13 - April 19) shifting 'b' up
% 280-480, 70 for 'a'
% 300-700, 50 for 'b'
% datapoint 2 is not a good inactivation
% datapoint 5 is a better inactivation
% 2. BK15BK14 experiment 2

% In the 4th experiment, 11 may be the last informative one


% PK39BK5 - note 'a' - [2000 2800]  something like 900-1200, 50

% PK16R47 - note 'e' [2800 3800] -- 240:440, b<60

% BK20BK45 - note 'a' [2800 2800] -- 220:420 b<60

    fvals=findwnoteJC('batch1002Anotes','a','','',0,[2000 3000],8500,1,'obs0',0);
        % 30ms backwards to note onset --- this measure of note onset hits too
        % early
    clear shifted
    for i=1:length(fvals)
        shifted(i,:)=fvals(i).datt;
    end
    
    
    [pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,2000,3000,[1 2 3],'obs0',1);
   %Exclude bad ones
   
    count=0;
    clear pitchFinal;
    for i=1:size(pitch,2)
        a=pitch(220:420,i);
        b=std(a);
        if b<60
            count=count+1;
            pitchFinal(:,count)=pitch(:,i);
        end
    end
    
    
    for i=4
        count=0;
        clear pitch;
        dept=0;
        pitch=BK20BK45_A(i).pitchcurves;
        clear BK20BK45_A(i).pitchcurves;
        for k=1:size(pitch,2)
            a=pitch(220:420,k);
            for j=2:length(a)
                b(i-1)=abs(a(i)-a(i-1));
            end
            if max(b)>25
                dept=1;
            end
            if dept==0
                count=count+1;
                BK20BK45_A(i).pitchcurves(:,count)=pitch(:,k);
            end
        end
    end

n=19;

BK20BK45_A(n).pitchcurves=pitch;
BK20BK45_A(n).rawdata=shifted;



    
[vals,trigs]=triglabel('batch1002Afiles','a',1,1,0,0);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

toffset=((toff/1000)*(32000)-512)/4+240;




BK20BK45_A(n).toffset=toffset;






B(1).baseline(n)=0;
B(1).acsf(n)=0;


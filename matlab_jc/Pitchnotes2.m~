% 1. BK15BK14 experiment 1 (April 13 - April 19) shifting 'b' up
    % 280-480, 70 for 'a'
    % 300-700, 50 for 'b'
    
   
% datapoint 2 is not a good inactivation
% datapoint 5 is a better inactivation

% 2. BK15BK14 experiment 2

% In the 4th experiment, 11 may be the last informative one

% BK15BK14_exp1 - [2000 3000] -- 300:650
% BK15BK14_exp2 - [3000 4100] -- 280:420
% BK15BK14_exp3 - [2200 2800] -- 280:580
% BK15BK14_exp4 - [2900 3900] -- 310:460

% PK39BK5 - note 'a' - [2000 2800]  -- 900:1200
% PK16R47 - note 'e' [2800 3800] -- 250:400

% BK20BK45 - note 'a' [2800 4000] -- 220:420
        %    note 'b' [2000 3000] -- 300:600
        
% PK20R49 - note 'a'  -- 260:410(?)
% PK32BK28_exp2 - note 'b' -- 400:700

% PU34 - note 'a' - long stack -- 750:1050

% BK50W18 - note 'a' - 500:750 - [2200 3000]
  % - for FrontDown - do 450:700


    fvals=findwnoteJC('batch1024Anotes','a','','',0,[2000 3000],8500,1,'obs0',0);
        % 30ms backwards to note onset --- this measure of note onset hits too
        % early
    clear shifted
    for i=1:length(fvals)
        shifted(i,:)=fvals(i).datt;
    end
    
    
    [pitch,avg]=jc_pitchmat1024(shifted,1024,1020,1,2000,3000,[1 2 3],'obs0',1);
   %Exclude bad ones
    
    
n=13;
% 
% B(n).pitchcurves=pitch;
% B(n).rawdata=shifted;



n=n+1;    
[vals,trigs]=triglabel('batch1024Afiles','a',1,1,0,0);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

toffset=((toff/1000)*(32000)-512)/4+240;

C(n).data=toffset;






B(1).baseline(n)=0;
B(1).acsf(n)=0;






B=BK15BK14_exp4;
A=jctester2(B,2000,3000);
A=jctester3(A,500,750);
for i=1:length(A)
    B(i).pitchcurves=0;
    B(i).selectedpitchcurves=A(i).selectedpitchcurves;
    B(i).rawpitchcurves=A(i).pitchcurves;
end



% For my birds
% 1. transfer fvals
% C(5).dat=fvals
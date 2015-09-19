for i=1:size(pitch4D,2)
    resid4D(:,i)=pitch4D(:,i)-mean(pitch4D')';
end
for i=1:size(pitch4Pd,2)
    resid4Pd(:,i)=pitch4Pd(:,i)-mean(pitch4Pd')';
end
for i=1:size(pitch4P,2)
    resid4P(:,i)=pitch4P(:,i)-mean(pitch4P')';
end
for i=1:size(pitch4UD,2)
    resid4UD(:,i)=pitch4UD(:,i)-mean(pitch4UD')';
end


AlldataBFlesion(4).pitchUDpost=pitch4P;
AlldataBFlesion(4).pitchUDpre=pitch4UD;
AlldataBFlesion(4).pitchDpre=pitch4D;
AlldataBFlesion(4).pitchDpost=pitch4Pd;
AlldataBFlesion(4).residDpost=resid4Pd;
AlldataBFlesion(4).residUDpost=resid4P;
AlldataBFlesion(4).residUDpre=resid4UD;
AlldataBFlesion(4).residDpre=resid4D;



pitch4P=jc_hilbertns(AlldataBFlesion(4).rawdataUDpost,1500,3000,44100);
pitch4UD=jc_hilbertns(AlldataBFlesion(4).rawdataUDpre,1500,3000,44100);
pitch4D=jc_hilbertns(AlldataBFlesion(4).rawdataDpre,1500,3000,44100);
pitch4Pd=jc_hilbertns(AlldataBFlesion(4).rawdataDpost,1500,3000,44100);
for i=1:8
    [xvals(i).dat,PSDBF.psdUDpre(i).dat]=jc_psd(AlldataBFlesion(i).residUDpre(AlldataBFlesion(i).onset:AlldataBFlesion(i).offset,:));
    [a,PSDBF.psdDpre(i).dat]=jc_psd(AlldataBFlesion(i).residDpre(AlldataBFlesion(i).onset:AlldataBFlesion(i).offset,:));
    [a,PSDBF.psdUDpost(i).dat]=jc_psd(AlldataBFlesion(i).residUDpost(AlldataBFlesion(i).onset:AlldataBFlesion(i).offset,:));
end

% normalized by bandwidth - PSD per ms
figure;hold on;
for i=1:8
    plot(xvals(i).dat,PSDBF.psdUDpre(i).dat./xvals(i).dat,'k')
    plot(xvals(i).dat,PSDBF.psdDpre(i).dat./xvals(i).dat,'b')
    plot(xvals(i).dat,PSDBF.psdUDpost(i).dat./xvals(i).dat,'r')
end
figure;hold on;
for i=1:8
    plot(xvals(i).dat,PSDBF.psdUDpre(i).dat./PSDBF.psdDpre(i).dat,'b')
    plot(xvals(i).dat,PSDBF.psdUDpre(i).dat./PSDBF.psdUDpost(i).dat,'r')
end

for i=1:size(B.pitchUDpost,2)
    B.residUDpost(:,i)=B.pitchUDpost(:,i)-mean(B.pitchUDpost')';
end
for i=1:size(pitch9INA,2)
    resid9INA(:,i)=pitch9INA(:,i)-mean(pitch9INA')';
end




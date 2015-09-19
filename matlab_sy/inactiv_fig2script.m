%inactiv_fig2script.m
%main panel is made with call to plothists, and call to plotline.

%raw data panel is made with calls to inactivplotoneday.m
%top pael pulled off illustrator file.
%to make acsf day 1, target note.
%ps is short for plotstruct.
figure
ax(1)=subplot(3,3,1);
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1:3];
ps.dyind=4
ps.notevl=1;
ps.mu='mu'
ps.ax=ax(1);
ps.plot='all';
inactivplot4(avls,graphvals,ps);
axis off;
box off;
ax(2)=subplot(3,3,4);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1];
ps.dyind=6
ps.notevl=1;
ps.mu='mu'
ps.ax=0;%now make muscimol run.
inactivplot4(avls,graphvals,ps)
axis off;
ax(3)=subplot(3,3,7);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1:3];
ps.dyind=6
ps.notevl=1;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)

ax(4)=subplot(3,3,2);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1:3];
ps.dyind=4
ps.notevl=2;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)
axis off;
ax(5)=subplot(3,3,5);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1];
ps.dyind=6
ps.notevl=2;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)
axis off;
ax(6)=subplot(3,3,8);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1:3];
ps.dyind=6
ps.notevl=2;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)

ax(7)=subplot(3,3,3);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1:3];
ps.dyind=13
ps.notevl=1;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)
axis off;

ax(8)=subplot(3,3,6);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1];
ps.dyind=15
ps.notevl=1;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)
axis off;
ax(9)=subplot(3,3,9);
%now make acsf, pitchshifted day 3
ps.ht=[0.7 0.7]
ps.line=0;
ps.tickht=.05;
ps.zrat=0;
ps.hist=[1:3];
ps.dyind=15
ps.notevl=1;
ps.mu='mu'
ps.ax=ax(3);%now make muscimol run.
inactivplot4(avls,graphvals,ps)

linkaxes(ax(1:3))
axes(ax(1))
axis([6300 7200 0 0.8])

linkaxes(ax(4:6))
axes(ax(4))
axis([1800 2700 0 0.8])

linkaxes(ax(7:9))
axes(ax(7))
axis([6600 7500 0 0.8])

for ii=[2 3  8 9]
    axes(ax(ii))
    plot([avls.initmean{1} avls.initmean{1}], [.0 0.7],'k--','Linewidth',1)
end  
for jj=[5]
    axes(ax(jj))
    plot([avls.initmean{2} avls.initmean{2}],[0 0.6],'k--','Linewidth',1)
end



figure


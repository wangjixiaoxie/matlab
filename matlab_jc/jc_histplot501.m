function jc_histplot501(Dvector,UDvector,Lvector)
figure; hold on
D=histc(Dvector,[-200:1:200]);
D=D/max(D);
UD=histc(UDvector, [-200:1:200]);
UD=UD/max(UD);
L=histc(Lvector, [-200:1:200]);
L=L/max(L);
stairs([-200:1:200],D,'r')
stairs([-200:1:200],UD,'b')
stairs([-200:1:200],L,'k')

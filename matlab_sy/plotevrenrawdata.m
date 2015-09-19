ind=11;
ax=subplot(411);


   strvl=[baspath fndat(ind).nm];
   cmd=['cd ' strvl]
   eval(cmd);
   mkfigs;

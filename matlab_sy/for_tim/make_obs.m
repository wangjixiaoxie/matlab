function make_obs(songfile,song)
% 
% function make_obs(songname,song)
% <songfile>:name of file to create
% <song>:song vector

     %create file name by parsing original file name
     songname = [songfile,'.cbin'];
     
     nmbrpnts = length(song);
     
     %create obs file
     obswrite(song,char(songname));
     
     %create rec file
     recname = [songfile,'.rec'];
     fid = fopen(char(recname),'a');
     dateinfo = ['File created ',date]; 
     beginfo = ['begin rec =      0 ms'];triginfo = ['trig time =      0 ms']; endinfo = ['rec end = ', num2str(round(length(song)/32)),' ms']; 
     freqinfo = ['ADFREQ =  32000']; sampinfo = ['Samples =  ', num2str(nmbrpnts)]; chaninfo = ['Chans =  1'];
     fprintf(fid,'%s\n',dateinfo,beginfo,triginfo,endinfo,freqinfo,sampinfo,chaninfo);
     fclose(fid);     

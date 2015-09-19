function []=MoveRght();
xl=xlim;
xlim(xl+diff(xl)*0.9);
return;

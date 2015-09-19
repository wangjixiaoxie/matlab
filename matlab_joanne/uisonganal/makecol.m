function  col_vect=makecol(vect);

%makes sure that argument is a column vector;


[nrows, ncols]=size(vect);


if min([nrows, ncols]) > 1
  disp('make_col: argument is not a vector');
else
  if ncols > nrows
    col_vect=vect';
  else
    col_vect=vect;
  end
end

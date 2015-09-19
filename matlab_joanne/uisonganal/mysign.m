function sgn = mysign(x)

sgn = sign(x);
sgn(find(sgn == 0)) = 1;


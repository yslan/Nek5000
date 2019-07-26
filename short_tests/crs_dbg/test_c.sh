make

np=1;./neklmpi c $np; echo "np = "$np;grep 'relres, rmin, rmax ' c.log.$np; grep relresv c.log.$np 
np=2;./neklmpi c $np; echo "np = "$np;grep 'relres, rmin, rmax ' c.log.$np; grep relresv c.log.$np
np=3;./neklmpi c $np; echo "np = "$np;grep 'relres, rmin, rmax ' c.log.$np; grep relresv c.log.$np
np=4;./neklmpi c $np; echo "np = "$np;grep 'relres, rmin, rmax ' c.log.$np; grep relresv c.log.$np

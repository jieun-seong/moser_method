# This is a shell archive.  Save it in a file, remove anything before
# this line, and then unpack it by entering "sh file".  Note, it may
# create directories; files and directories will be owned by you and
# have default permissions.
#
# This archive contains:
#
#	cohom.m
#	der.m
#	initial
#	initialize.m
#	myshift.m
#	newton.m
#	sobolev.m
#	sobolev-high.m
#
echo x - cohom.m
sed 's/^X//' >cohom.m << 'END-of-cohom.m'
Xfunction cohomology = cohom(f)
X	 global cohomarray;; 
X	 cohomology = ifft(cohomarray.*fft(f));
Xendfunction
END-of-cohom.m
echo x - der.m
sed 's/^X//' >der.m << 'END-of-der.m'
Xfunction derivative = der(f); 
X	 global derarray; 
X	 derivative = ifft(derarray.*fft(f));
Xendfunction
X
END-of-der.m
echo x - initial
sed 's/^X//' >initial << 'END-of-initial'
X
X#{
X# This program implements a method in J. Moser 66 lectures 
X# to compute the conjugacy of a circle map to a rotation. 
X#
X# It is a program designed to illustrate the algorithm.  It is 
X# not a production quality program.
X#      It does not include checks of correctness of the programs
X#      The calculations proceed without checking that the results 
X#          are correct. It would suffice to estimate sobolev norms. 
X#      It uses complex Fourier transforms. Using real, one can reduce 
X#         the space and computation time by a factor 2, 
X#       We do not attempt to make continuation.
X#	We do not include the ability of increasing the number of terms 
X#	   dynamically. Changing N requires to restart and reinitialize.  
X#      The program wil develop instabilities if N  is bigger that 1024 or so. 
X#
X#      We have used mainly vector operations in Octave which 
X#      are comparable or faster than in C. So, the program is in 
X#      the ballpark of C programs.  It has been tested slightly in Matlab.
X#  
X# We consider f(x) = x + lambda*sin(x) + drift
X# We want to find h(x) = x + hath(x) (with hath periodic of period 2 pi) 
X# and drift  in such a way  that 
X# f (h(x)) + drift = h(x + twopiomega) 
X# 
X# Note we obtain a quadratically convergent method, but no matrix NxN is 
X# stored. Storage required is a multiple of N. 
X#}
X
X# N is the number of Fourier modes. 
XN = 1024; 
Xomega = (sqrt(5) -1)/2;
Xglobal twopiomega = 2*pi*omega
Xglobal lambda = 0.1;
X
Xglobal xarray;
Xglobal farray;  
Xglobal derarray; 
Xglobal shiftarray; 
Xglobal cohomarray; 
X
Xinitialize(N)
X
X# This is an inialization of the variables. 
Xhath = linspace(0,0, N); 
Xglobal drift  = twopiomega; 
X
X
X
X
END-of-initial
echo x - initialize.m
sed 's/^X//' >initialize.m << 'END-of-initialize.m'
Xfunction dummy = initialize(N) 
X# We set an array which contains x 
Xglobal xarray; 
Xglobal farray; 
Xglobal farranyhmask; 
Xglobal derarray; 
Xglobal shiftarray; 
Xglobal cohomarray; 
Xglobal twopiomega; 
X
X
X
Xxarray = linspace(0,2*pi*(1 -1/N), N);
X#The following lines set an array that tells which Fourier
X#coefficient appears at a place in FFT (read the manual!!) 
X# 
X# The derivative is obtained by multiplying that by I 
X# The shift is an exponential.
X
X
Xtempa = linspace(0,N/2 -1, N/2); 
Xtempb = linspace(N/2, N-1, N/2) - N; 
Xfarray = [tempa,tempb]; 
X
Xderarray = farray*I; 
X
Xtempa = linspace(0,0,N/4); 
Xtempb = linspace(1,1, N - 2*(N/4)); 
Xfarrayhmask = [tempa, tempb, tempa]; 
X
Xshiftarray = exp( twopiomega*derarray);
X
X
X# We want to get roughly cohomarray = 1./(shiftarray -1) but we 
X# have to take care that the coefficient 1 is 0.  We want to do 
X# this using only vector operations and not write loops. 
Xtemp = 1-shiftarray; 
Xtemp(1) = 1;  
Xtemp2 = 1./temp;
Xtemp2(1) = 0; 
Xcohomarray = temp2;
X
X
X
Xendfunction
X
X
X
END-of-initialize.m
echo x - myshift.m
sed 's/^X//' >myshift.m << 'END-of-myshift.m'
Xfunction shifted = myshift(f)
X	 global shiftarray; 
X	 shifted =  ifft(shiftarray.*fft(f));
Xendfunction
X
END-of-myshift.m
echo x - newton.m
sed 's/^X//' >newton.m << 'END-of-newton.m'
Xfunction newhath  = newton(hath)
Xglobal lambda;
Xglobal drift; 
Xglobal xarray;
Xglobal twopiomega;
X# This algorithm is explained in J. Moser's Pisa Lectures in 66.
X# It uses cancellations based on what J. Moser calls "group structure". 
X# Note that we do not use composition. 
X#
X#
XerrorI = hath + lambda*sin(xarray + hath) -  myshift(hath)+drift+twopiomega;
Xnorm(errorI)
Xderh  = der(hath) +1;
Xmyfactor = myshift(derh); 
Xdelta = - (sum(errorI./myfactor))/(sum(1.0./myfactor));
Xw =  cohom( (-errorI + delta)./myfactor); 
Xdrift = drift +delta; 
Xnewhath=derh.*w + hath;
Xendfunction
END-of-newton.m
echo x - sobolev.m
sed 's/^X//' >sobolev.m << 'END-of-sobolev.m'
Xfunction sobnorm =sobolev(s, fhat); 
X	 global derarray; 
X        sobnorm = norm( (derarray.^s).*fhat); 
Xendfunction
X
END-of-sobolev.m
echo x - sobolev-high.m
sed 's/^X//' >sobolev-high.m << 'END-of-sobolev-high.m'
Xfunction sobnormh =sobolevh(s, fhat); 
X  global derarray; 
X  global farrayhmask
X
X
X  sobnormh = norm( (derarray..^s)..*fhat.*farrayhmask); 
Xendfunction
X
END-of-sobolev-high.m
exit


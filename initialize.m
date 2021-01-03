function dummy = initialize(N) 
# We set an array which contains x 
global xarray; 
global farray; 
global farrayhmask; 
global derarray; 
global shiftarray; 
global cohomarray; 
global twopiomega; 



xarray = linspace(0,2*pi*(1 -1/N), N);
#The following lines set an array that tells which Fourier
#coefficient appears at a place in FFT (read the manual!!) 
# 
# The derivative is obtained by multiplying that by I 
# The shift is an exponential.


tempa = linspace(0,N/2 -1, N/2); 
tempb = linspace(N/2, N-1, N/2) - N; 
farray = [tempa,tempb]; 

derarray = farray*I; 

tempa = linspace(0,0,N/4); 
tempb = linspace(1,1, N - 2*(N/4)); 
farrayhmask = [tempa, tempb, tempa]; 

shiftarray = exp( twopiomega*derarray);


# We want to get roughly cohomarray = 1./(shiftarray -1) but we 
# have to take care that the coefficient 1 is 0.  We want to do 
# this using only vector operations and not write loops. 
temp = 1-shiftarray; 
temp(1) = 1;  
temp2 = 1./temp;
temp2(1) = 0; 
cohomarray = temp2;



endfunction




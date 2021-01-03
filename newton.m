function newhath  = newton(hath)
global lambda;
global drift; 
global xarray;
global twopiomega;
# This algorithm is explained in J. Moser's Pisa Lectures in 66.
# It uses cancellations based on what J. Moser calls "group structure". 
# Note that we do not use composition. 
#
#
errorI = hath + lambda*sin(xarray + hath) -  myshift(hath)+drift+twopiomega;
norm(errorI)
derh  = der(hath) +1;
myfactor = myshift(derh); 
delta = - (sum(errorI./myfactor))/(sum(1.0./myfactor));
w =  cohom( (-errorI + delta)./myfactor); 
drift = drift +delta; 
newhath=derh.*w + hath;
endfunction

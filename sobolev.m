function sobnorm =sobolev(s, fhat); 
	 global derarray; 
        sobnorm = norm( (derarray.^s).*fhat); 
endfunction


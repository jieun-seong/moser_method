function cohomology = cohom(f)
	 global cohomarray;; 
	 cohomology = ifft(cohomarray.*fft(f));
endfunction

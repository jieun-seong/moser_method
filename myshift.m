function shifted = myshift(f)
	 global shiftarray; 
	 shifted =  ifft(shiftarray.*fft(f));
endfunction


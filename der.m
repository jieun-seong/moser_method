function derivative = der(f); 
	 global derarray; 
	 derivative = ifft(derarray.*fft(f));
endfunction


function sobnormh =sobolevh(s, fhat); 
  global derarray; 
  global farrayhmask


  sobnormh = norm( (derarray..^s)..*fhat.*farrayhmask); 
endfunction


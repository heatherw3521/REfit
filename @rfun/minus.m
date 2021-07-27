 function h = minus(r, g)
 % subtract g from r and represent the result as an rfun. 
 % One of r or g must be an rfun. The other can be a scalar or an efun.  
 %
 % minus(r, g) is called for synax 'r-g'. 
 %
 %%
 % See also rfun/plus, rfun/compress. 
 g = uminus(g); 
 h = plus(r, g);
 end
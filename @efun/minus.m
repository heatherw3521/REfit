 function h = minus(s, g)
 % subtract g from s and represent the result as an efun. 
 % One of s or g must be an efun. The other can be a scalar or an efun.  
 %
 % minus(s, g) is called for synax 's-g'. 
 %
 %%
 % See also efun/plus, efun/compress. 
 g = uminus(g); 
 h = plus(s, g);
 end
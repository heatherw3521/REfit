function p = poles(s, varargin)
% poles(s) returns the poles of the trigonometric rational 
% r(x) = ift(s). 
%
% poles(s, 'zt') returns the poles of the rational r(z), where
% z = exp(2*pi*1i*(x-s.domain(1))/diff(s.domain)). 
%
%%
% See also: efun/roots. 

if isempty(s)
    p = [];
    return
end
a = s.domain(1); 
b = s.domain(2); 
per = b-a; 

if strcmp(s.space, 'value')
    error('efun:getpoles:cannot get poles for an efun in time/signal space')
else   
    p = exp((s.exp)); 
  
    p = [p; 1./conj(p)]; 
    if isempty(varargin) %poles wrt to x, not z. 
    % adjust for domain:
        p = -log(p)/2/pi/1i*per;
            %adjust to interval [a,b)
        p = a + mod(real(p), per) + imag(p)*1i;
    %adjust to domain: 
    %rp = (1-a)/b*real(p)+a; 
    %rp = real(p) + a; 
    %p = rp + imag(p)*1i;  
    %if 
        %p = exp(-2*pi*1i*(p-a)/per); 
        %p = p(abs(p)<1); %only return poles inside disk.
    end
end

end


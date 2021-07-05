function p = poles(s, varargin)
% poles(s) returns the poles of the trigonometric rational r(x) = ift(s). 
%
% poles(s, 'zt') returns the poles of the rational r(z), where r(z) = r(x)
% on the unit circle and z = exp(2*pi*1i*x). 
%
% See also: roots. 

if isempty(s)
    p = [];
    return
end

if strcmp(s.space, 'value')
    error('efun:getpoles:cannot get poles for an efun in time/signal space')
else   
    p = exp((s.exp)); %these are the poles wrt z = e^(-i*2*pi theta)
    % adjust for domain: 
    a = s.domain(1); 
    b = s.domain(2); 
    p = [p; 1./conj(p)];  
    p = -log(p)/2/pi/1i;
        %adjust to interval [0,1)
    p = mod(real(p), 1) + imag(p)*1i;
    %adjust to domain: 
    rp = (1-a)/b*real(p)+a; 
    p = rp + imag(p)*1i;  
    if ~isempty(varargin) && strcmpi(varargin{:}, 'zt')
        p = exp(-2*pi*1i*p); 
        p = p(abs(p)<1); %only return poles inside disk.
    end
end

end


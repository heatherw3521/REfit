%EXAMPLE SCRIPT 1: Introduction to REfit. 

% The REfit package allows one to construct and compute with
% trigonometric rational functions fit to periodic or compactly supported 
% data in time/signal space, and correspondingly, fits by sums of decaying 
% exponentials in Fourier space. See [1,4] for more information. 
%
% Trigonometric rationals are encoded as rfun objects, and sums of exponentials 
% in Fourier space are encoded as efun objects. The syntax for computing with
% efuns and rfuns is largely inspired by Chebfun [5]. Many standard MATLAB 
% commands are overloaded to make computing with these objects easier. 
%
% This example file demonstrates how to construct efuns and rfuns, and
% how to use some of the basic commands associated with them.
% 

%%
% We start with a function handle and a grid for sampling the function. 
clear all
n = 800;  
f = @(x) (abs(sin(2*pi*x)).^3+ exp(sin(2*pi*x)))/4; 
x = linspace(0,1,2*n+2); x = x(1:end-1).';

%%
% We can construct an approximation to the function directly using
% rfun:

r = rfun(f(x), x) 
plot(r); shg

%%
% The rfun object r encodes a type (45, 46) barycentric trigonometric rational
% interpolant [2,4] to the sample, along with a constant term. Without the
% constant term, r is approximately a mean zero function. Stored in r are 
% the components making up the interpolant.  

[r.nodes,...     % Barycentric nodes (interpolating points)
r.vals,...      % value of interpolant at the nodes.
r.weights]   % Barycentric weights. 
r.const     % constant = mean(r) on [0, 1).
%%
% By default, the domain of r is given by [0, 1).

r.domain

% This can be changed to [a, b] in the call to rfun using the syntax 
% rfun(r, 'domain', [a, b]). It is assumed that r is periodic with 
% period (b-a). 
%%
% The rfun constructor tries to resolve the sample with a default tolerance 
% parameter of 1e-9. The tolerance can be adjusted by typing 
% rfun(f(x), x, 'tol', tol), where tol is the desired tolerance. 

clf
semilogy(x,abs(r(x)-f(x)), 'linewidth', 1.5)
title('absolute error')
shg

%%
% r has 46 poles. In infinte precision, they occur in conjugate pairs.
% We plot them below. 
clf
plot(poles(r), '.k', 'markersize', 20)
hold on
plot([0 1], [0 0], '--k')
hold off
title('Poles of r in complex plane')
shg

%%
% For a quick visual of poles and their locations relative to the
% features of r, we use the poleplot command. This displays the 
% poles of r in the complex plane, superimposed over a plot of r on [0, 1).  

clf; 
poleplot(r)
title('poles of r with r(x) superimposed on complex plane')

% Notice that the poles cluster toward the singularities of the function f. 
% Rfun is based on the AAA algorithm. For more information about AAA and 
% pole-clustering, see [3]. 
%%
% The poles above are given with respect to r(x), where x is on [0, 1]. One
% may also be interested in the poles of the rational function r(z), where
% z = exp(2*pi*1i*x). These are accessed with the 'zt' (z-transform) flag:

pz = poles(r, 'zt'); 
clf
poleplot(r, 'zt')
title('poles of r(z) in unit disk')
shg
% Poleplot only plots the poles of r(z) that occur inside the unit
% disk.

%%
% Several MATLAB and Chebfun-like commands [5] are available. For example,
% we can take the indefinite integral of r. By default, it is returned as a 
% function handle.
  
R = cumsum(r); 
R(.35)
tru = 0.23455361203219054
abs(R(.35)-tru)

% We can also specify that R should be returned as an rfun, but
% since rfuns are periodic, this only makes sense if r is a mean zero 
% function. The result is the periodic extension of cumsum(r) on [0, 1). 
% Here, we plot such a result on [0, 2). 
r = r-mean(r); 
R = cumsum(r, 'rfun'); 
plot(linspace(0, 2, 500).', R(linspace(0,2, 500).'), 'r')

%%
% We can differentiate R to get r back. We specify that we'd like an
% rfun representation of the derivative, rather than just a function handle. 

r2 = diff(R, 'rfun'); 
clf
semilogy(abs(r2(x)-r(x)))

%% 
% Here, we plot the roots and local minima/maxima of r: 

rt= roots(r); 
mins = min(r); 
mxs = max(r); 

clf
plot(rt, r(rt), '.k', 'markersize', 20)
hold on
plot(mins, r(mins), '.r', 'markersize', 20)
plot(mxs, r(mxs), '.b', 'markersize', 20)
plot(x, r(x), 'linewidth', 2);
legend('roots', 'minima', 'maxima') 
title('peaks, valleys, and roots of r')
hold off
shg

%%
% Another way to fit a trigonometric rational to a given sample
% is by working in Fourier space. The Fourier transform of a type (k-1, k)
% trigonometric rational function is a sum of decaying complex
% exponentials [1,4]. The efun class in REfit constructs sums of exponentials
% that approximately interpolate the Fourier coefficients associated with a
% signal. 

%% 
% To illustrate the idea, consider a cubic spline on [0, 1). 
% We retrieve a function handle and plot the spline in signal space
% below, making use of a gallery function where the spline information is
% stored for convenience. 

[~, fs] = gallery_efun('spline'); 
clf
plot(x, fs(x), 'linewidth', 2)
title('a cubic spline')

%%
% We can call efun using a sample of fs. The constructor uses the fft
% to find the Fourier coefficients associated with a bandlimited
% approximation to fs, where the bandlimit = floor((M-1)/2), M =
% length(fs(x)). 

S = efun(fs(x))

%%
% The efun object S represents an exponential sum of the form 
% S(k) = sum_{j = 1}^{} w_j exp(lambda_j k), where the real part of
% lambda is negative. For k = 0, 1, ..., the value S(k) is approximately 
% equal to the kth Fourier coefficient of the cubic spline. 
%%
S(0) % the Fourier coefficient of fs at the zeroth mode is 
     % given by evaluating s at zero.
     
mean(fs(x)) % the mean value of fs (also the value of the zeroth Fourier coeff). 

%%
% The accuracy of S as a representation of the cubic spline depends
% on the number of samples [4]. We store this information as a 'resolution
% parameter: 

N = S.res
%%
% S serves as a representation of a bandlimited projection of the 
% cubic spline. Here we plot the error in Fourier space within the
% bandlimit.

coeffs = sample2coeffs(fs(x)); %Fourier coeffs of fs, via fft. 
 
clf
semilogy(-N:N, abs(coeffs-S(-N:N)))
title('Absolute error: S vs. fft')

%%
% efun can also be constructed using Fourier data. We include the
% modes and a special 'coeffs' flag.
S = efun(coeffs, -N:N, 'coeffs')

%%
% We can plot the rational given by the inverse Fourier transform of S
% using the 'plot' command:
clf
plot(S)
title('efun representation of cubic spline')

%%
% We can quickly plot Fourier coefficients as well
clf
plotcoeffs(S)

% (This command also works for rfuns.) 
%%
% As with rfun, by default, the signal domain is assumed to be [0, 1) and
% the tolerance is set to 1e-9. These parameters can be changed using the 
% syntax efun(...,'tol', tol, 'domain', [a,b])
%
% The poles of the rational associated with the inverse Fourier transform
% of s are also easily accessible. One can see that they cluster around
% the knots of the cubic spline:

clf
plot(poles(S), '.r', 'markersize', 20)
hold on
plot([1, 2, 3, 4, 5]/6, zeros(5,1), '.k', 'markersize', 20)
plot([0, 1], [0, 0], '--k')
legend('poles', 'knot locations')
title('poles of the inverse Fourier transform of s')

%%
% The efun format is closely connected to the pole-residue form of the inverse
% Fourier transform of S. We can access a handle for the pole-residue
% form using the ift function: 

h = ift(S, 'polres'); 
clf
plot(x, h(x), 'linewidth', 2)
title('Evaluation of an efun in signal space')

%%
% we can also evaluate the inverse Fourier transform of S on its domain using
% the following syntax, where the 'values' flag indicates that we want to 
% evaluate S in value/signal space, rather than Fourier space:
hold on
vals = S(x, 'values'); 
plot(x, vals, 'k', 'linewidth', 2)
hold off
%%
% Sometimes it is desirable to transform from an efun representation to
% an rfun representation of the rational function encoded in h. 
% For this, we use the ift command without any flags: 

s = ift(S) %the rfun s. 
clf
poleplot(s)

%%
% In the same way, we can convert rfuns to efuns with the "Fourier
% transform" command ft:
clf
R = ft(r) %constructing an efun R
poleplot(R)
%%
% For information on how rfuns and efuns can be used together to overcome
% approximation challenges, and/or to explore the advantages/disadvantages 
% of each representation, see [4], as well as the example at the end of 
% this script.

%%
% Now we have efuns S and R, and correspondingly, rfuns s and r. Objects of
% the same type can be added, multiplied, or convolved. 

RS = (R-.3*S).*(S) %combining efuns with sums and products
clf
plot(R-.3*S)
hold on
plot(S+.2)
plot(RS)
legend('R-.3*S', 'S', 'convolution')
title('products of efuns = convolutions of signals they represent')
%%
% The product of efuns in Fourier space is equivalent to the convolution of 
% rfuns in signal space. 
rs = convolve(r-.3*s, s); 
max(abs(RS(x, 'values')-rs(x)))

%%
% It is often the case that when one form of approximation struggles, the other
% method can succeed. In the next example, rfun struggles to return a
% stable representation for the product of two rfuns. 

rps = (r-.3*s).*(s)
clf
plot(r-.3*s)
hold on
plot(s)
plot(rps)
legend('r-.3s', 's', 'product')
title('product of r-.3*s and s')
%%
% The rational approximation looks fine to the naked eye, but it contains so--called
% 'spurious poles' or "Froissart doublets" [2,4] that occur on the real line. 

find(abs(imag(poles(rps))) < 1e-11)
%%
% The presence of these poles can lead to evaluation errors; It is better
% if we can construct a representation that avoids them. We try using
% efuns instead. 

RPS = convolve(R-.3*S, S) %times in signal space = convolve in Fourier space
clf
plot(R-.3*S)
hold on
plot(S)
plot(RPS)
legend('r-.3s', 's', 'product')
title('product of signals assoc. with R-.3*S and S')

%%
% The efun poles are certain to occur in conjugate pairs and will always have 
% nonzero imaginary parts [1,4]. If the rfun format is still desirable, 
% we can use the ift to convert to an rfun. This approach involves a very 
% different approximation process. As we verify below, the resulting rfun 
% does not have any real-valued poles. 

rps = rfun(RPS)
clf
semilogy(x,abs(RPS(x, 'values')-rps(x)), 'linewidth', 1.5)
title('absolute error in efun and rfun representations of the product')

% check for real-valued poles:
find(abs(imag(poles(rps))) < 1e-11)

%%
% There are many more functions to explore (and more to be developed!).
% Please email heatherw3521@gmail.com if you have ideas for functions, 
% examples, or are curious about how we can apply efuns and/or rfuns in
% your work.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% REFERENCES%%%%%%%%%%%%%%%%%%

% [1] Beylkin, G., and Monzón, L., "Nonlinear inversion of a band-limited 
% Fourier transform", Appl. and Comp. Harm. Anal., (27) 351--366, (2009).
%
% [2] Nakatsukasa, Y., Sète, O., and Trefethen, L.N, "The AAA algorithm
% for rational approximation", SIAM J. Sci. Comp., 40(3), A1494–-A1522,
% (2018).
%
% [3] Trefethen, L.N., Nakatsukasa, Y., and Weideman, J.C., "Exponential 
% node clustering at singularities for rational approximation, quadrature, 
% and PDEs", Num. Math., (147) 227-–254, (2021).
%
% [4] Wilber, H., Damle, A., and Townsend, A., "Data-driven algorithms for 
% signal processing with rational functions", arXiv:2105.07324, (2021).
%
% [5] www.chebfun.org


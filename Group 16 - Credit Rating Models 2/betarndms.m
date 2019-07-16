function r = betarndms(mu, sigma, m, n)
%BETARNDMS Random arrays from beta distribution given the mean and standard
%deviation.
%   R = BETARNDMS(MU,SIGMA) returns a single number drawn from a beta
%   distribution with mean MU and standard deviation SIGMA.
%
%   R = BETARNDMS(MU,SIGMA,M,N) or R = BETARNDMS(MU,SIGMA,M) returns an
%   M-by-N or M-by-M array, respectively.
%
%   See also BETARND

%  Closed-form solution for the beta parameters a and b given mu and sigma:
a = (mu.^2 - mu.^3 - mu.*sigma.^2) ./ sigma.^2;
b = (mu - 2*mu^2 + mu.^3 - sigma.^2 + mu*sigma.^2) ./ sigma.^2;

if nargin == 2
    m = 1; n = 1;
elseif nargin == 3
    n = m;
end

r = betarnd(a, b, m, n);
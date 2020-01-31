function part = mutation(part,t)
% The function "mutation" mutates the particles according to a discrete 
% stochastic nonlinear dynamical systems. 
% 
% x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
% 
% Inputs: 
%   - part : particles
%   - t : time 
% 
% Output: 
%   - part : mutated particles 
% 
% Implemented example: 
%   x(t+1) = x(t)/2 + 25*x(t)/(1+x(t)^2) + 8*cos(1.2t) + w(t)
%   w(t) ~ N(0,1)
%
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

n_part = size(part,2);
part = part/2 + 25*part./(1+part.^2) + 8*cos(1.2*(t)) + randn(1,n_part)*sqrt(1);
end
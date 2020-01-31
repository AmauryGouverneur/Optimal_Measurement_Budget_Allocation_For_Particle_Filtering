function x = model(T)
% The function "model" simulates a realisation of a state vector according 
% to a discrete stochastic nonlinear dynamical systems. 
% 
% x(0) ~ F
% x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
%
% where w_t is the process noise. 
%
% Input: 
%   - T : length of the time interval  
% 
% Output: 
%   - x : state vector 
%
% Implemented example:
%   x(t+1) = x(t)/2 + 25*x(t)/(1+x(t)^2) + 8*cos(1.2t) + w(t)
%   x(0) ~ N(0,5^2)
%   w(t) ~ N(0,1)
% 
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

x = zeros(1,T+1);
x(0+1) = randn()*5;

for t = 0:T-1
    index_t = t+1;
    w_t =  randn()*sqrt(1);
    x(index_t+1) = x(index_t)/2 + 25*x(index_t)./(1+x(index_t).^2) + 8*cos(1.2*(t))+w_t;
end

end
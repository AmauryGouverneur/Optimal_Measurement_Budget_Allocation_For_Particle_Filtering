function y_j = measurements(x_j,T)
% The function "measurements" simulates a realisation of a measurement 
% vector of x_j according  to a discrete stochastic nonlinear dynamical
% systems. 
% 
% y(t) = g(x(t),v(t)) 
%
% Inputs: 
%   - x_j : state vector
%   - T : length of the time interval
% 
% Output: 
%   - y_j : observation vector 
%
% Implemented example: 
%   y(t) = (x(t)^2)/2 + v(t) 
%   v(t) ~ N(0,(sin(0.25*(0:T))+2)^2)
% 
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel 

v_j = randn(1,T+1).*(sin(0.25*(0:T))+2);
y_j = x_j.^2/20 + v_j ; 

end

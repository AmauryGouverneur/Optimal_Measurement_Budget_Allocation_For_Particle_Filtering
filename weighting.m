function w = weighting(y_j,part,t)
% The function "weighting" return the probability of observing y_j(t) given
% part is the "true" state vector at time t
%
% y(t) | x_(t) ~ p(y(t) | x(t)) dy(t)  
%
% Input: 
%   - y_j : measurement vector
%   - part : particles 
%   - t : time 
% 
% Outputs : 
%   - w : weigth vector
% 
% Implemented example: 
%   z(t) = x(t) 
% 
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel


index_t = t+1;
w = normpdf(y_j(index_t),(part.^2)./20,(sin(0.25*(t))+2))';

end
function [tau_j,part] = particle_filter(y_j,measurement_times,n_part,T) 
% The function "particle_filter" computes an estimate of the state vector
% x_j. This is performed using a SISR particle filter using measurement y_j
% given at a given measurements times. 
% 
% A discrete stochastic nonlinear dynamical system is modelled:
%   x(t+1) = f(x(t),w(t)) for t = 0,...,T ? 1,
%   y(t) = g(x(t),v(t)) for t in a given measurment set
%   z(t) = h(x(t)) for t = 0,...,T 
%   x(0) ~ F 
%   w(t) and v(t) are the process and measurement noise respectively witn
%   known probability density functions
%
% The following functions are implemented according to the discrete 
% stochastic nonlinear dynamical system : 
%   - initialization(n_part) : initializes the particles according to F 
%   - weighting(y_j,part,t) : return p(y_j|part) at time t
%   - mutation(part,t): mutates the particles according to the dynamical 
%   system
% 
% Inputs: 
%   - y_j : observation vector
%   - measurement_times : binary vector of size T+1, 1 indicating a 
%   measurment time, otherwise 0 
%   - n_part : number of particles
%   - T : length of the time interval
% 
% Outputs: 
%   - tau_j : estimate of x using the particle filter
%   - part : particles at the last time step 
%
% Date : 30/01/20
% Author : Amaury Gouverneur & Antoine Aspeel

part = initialization(n_part) ;
size_part = size(part,1);
tau_j = zeros(size_part,T+1);

if measurement_times(0+1)
        w = weighting(y_j,part,0); % weighting
        tau_j(:,1) = part*w/sum(w); % estimation
        ind = randsample(n_part,n_part,true,w);
        part = part(:,ind);
else 
    tau_j(:,1) = sum(part,2)/n_part;
end

for t = 1:T
    %mutation 
    part = mutation(part,t-1);
    index_t = t+1; 
    if measurement_times(t+1)
        w = weighting(y_j,part,t);
        tau_j(:,index_t) = part*w/sum(w); 
        %selection
        if norm(w)==0 
            tau_j = NaN; 
            break;
        else
            ind = randsample(n_part,n_part,true,w);
            part = part(:,ind);
        end
    else
        tau_j(:,index_t) = sum(part,2)/n_part;
    end 
end

end
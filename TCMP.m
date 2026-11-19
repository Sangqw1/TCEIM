function [a0,b0,aM,K] = TCMP(u,q_prev,gap,pGG,pTT,nu)
% TCMP  Temporally Correlated Mixture Prior.
% Inputs: normalized innovations u, previous Gaussian probabilities q_prev,
% elapsed steps gap, transition probabilities pGG/pTT, Student-t dof nu.
% Outputs: Beta shapes a0/b0, Markov prediction aM, adaptive gain K.

%% 1. Align channels and check the inputs
u = u(:);
q_prev = q_prev(:);
gap = gap(:);
assert(~isempty(u) && numel(u)==numel(q_prev) && numel(u)==numel(gap));
assert(all(isfinite(u)));
assert(all(isnan(q_prev) | (isfinite(q_prev) & q_prev>=0 & q_prev<=1)));
assert(isscalar(pGG) && isfinite(pGG) && pGG>0 && pGG<1);
assert(isscalar(pTT) && isfinite(pTT) && pTT>0 && pTT<1);
assert(isscalar(nu) && isfinite(nu) && nu>0);
known = isfinite(q_prev);
assert(all(isfinite(gap(known)) & gap(known)>=1 & gap(known)==round(gap(known))));

%% 2. Predict the Gaussian probability using the Markov model
% Stationary probability and temporal decay factor:
% piG = (1-pTT)/(2-pGG-pTT), phi = pGG+pTT-1.
piG = (1-pTT)/(2-pGG-pTT);
phi = pGG+pTT-1;

% New channels start from piG; known channels propagate their last posterior.
% aM = piG + phi^gap * (q_prev-piG).
aM = piG*ones(size(u));
aM(known) = piG + phi.^gap(known).*(q_prev(known)-piG);
aM = min(max(aM,eps),1-eps);

%% 3. Calculate the current Gaussian/Student-t evidence
% Both component densities use the same unit scale.
z = u.^2;
logG = -0.5*log(2*pi)-0.5*z;
logT = gammaln((nu+1)/2)-gammaln(nu/2)-0.5*log(nu*pi) ...
    -0.5*(nu+1)*log1p(z/nu);
logBF = logG-logT;
logBF(isinf(z)) = -Inf;

% c = BF/(1+BF), evaluated as a numerically stable sigmoid.
c = exp(-max(0,-logBF))./(1+exp(-abs(logBF)));

%% 4. Compute the adaptive evidence weight
% Bernoulli variance proxies; these are not calibrated KF error variances.
P = aM.*(1-aM);
R = c.*(1-c);
K = P./(P+R);

%% 5. Form the unit-concentration Beta initialization
% a0 = aM + K*(c-aM), b0 = 1-a0.
fused_probability = aM + K.*(c-aM);
a0 = min(max(fused_probability,eps),1-eps);
b0 = 1-a0;
end

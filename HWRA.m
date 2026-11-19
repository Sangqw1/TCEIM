function [p_bar,w] = HWRA(q,keep,p_prior)
% HWRA  Heavy-Tail-Weighted Risk Allocation.
% Inputs: converged Gaussian probabilities q, retained-channel matrix keep,
% and original hypothesis probabilities p_prior (H0 is the first row).
% Outputs: reallocated probabilities p_bar and conditional fault weights w.

%% 1. Align channels and check the hypothesis representation
q = q(:);
p_prior = p_prior(:);
assert(~isempty(q) && ~isempty(p_prior));
assert(size(keep,2)==numel(q) && size(keep,1)==numel(p_prior));
assert(all(isfinite(q) & q>=0 & q<=1));
assert(all(keep(:)==0 | keep(:)==1) && all(keep(1,:)==1));
assert(all(isfinite(p_prior) & p_prior>=0 & p_prior<=1));
assert(sum(p_prior)<=1+1e-10);

%% 2. Extract heavy-tail probabilities and the fixed monitored fault mass
rho = 1-q;
mass = sum(p_prior(2:end));              % P(M), excluding H0
p_bar = p_prior;
w = zeros(size(p_prior));
if mass==0
    return;
end

% Only positive-prior fault hypotheses participate in redistribution.
rows = find(p_prior>0);
rows(rows==1) = [];

%% 3. Score each distinct physical fault set
% s_j = product(rho_i over excluded channels) * product(q_i over retained).
% Identical exclusion patterns are scored once, avoiding duplicate counting.
[patterns,~,group] = unique(keep(rows,:),'rows');
logQ = log(max(q,realmin));
logRho = log(max(rho,realmin));
logS = patterns*logQ + (1-patterns)*logRho;

%% 4. Normalize within the monitored fault event
% Subtracting the largest log score prevents numerical overflow.
weights = exp(logS-max(logS));
weights = weights/sum(weights);

% Split each physical-set weight equally among its duplicate labels.
counts = accumarray(group,1);
w(rows) = weights(group)./counts(group);

%% 5. Restore the original total probability mass
% p_bar_j = P(M)*w_j. H0 and zero-prior fault rows retain their values.
p_bar(2:end) = mass*w(2:end);
assert(p_bar(1)==p_prior(1));
assert(abs(sum(p_bar(2:end))-mass)<=1e-12*mass);
end

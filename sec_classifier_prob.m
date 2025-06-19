function [C_sp, C_ss] = sec_classifier_prob(i_max,rp,rs,ks,p_obs)


%size_prob(i,n) = probability that a cluster with i primaries has overall
%size n
size_prob = zeros(i_max);
for i = 1:(i_max-1)
    for n = (i+1):i_max
        logl_ni = log(i) - log(n) + gammaln(ks*n+n-i) ...
            - gammaln(n-i+1) - gammaln(ks*n) ...
            + (n-i)*log(rs/ks) - (ks*n+n-i)*log(1+rs/ks);
        size_prob(i,n) = exp(logl_ni);
    end
end

%p_pdf = probability distribution for primary cases in a cluster
i_arr = 1:i_max;
p_pdf = (1 - 1/rp).^(i_arr-1)/rp;

%c_pdf(i,n) = probability that a cluster has i primaries and total size n;
c_pdf = size_prob.*repmat(p_pdf,i_max,1)';

%r_pdf = probability of a secondary case being the jth case in a cluster
r_pdf = zeros (1,i_max);
for j = 2:i_max
    sub_matrix = c_pdf(1:(j-1),j:i_max);
    r_pdf(j) = sum(sub_matrix(:));
end
%r_pdf needs normalization
r_pdf = r_pdf*(1-rs)/rp/rs;

if sum(r_pdf) < .98
    % 'sum r_pdf is insufficient'
    % [i_max rp rs ks p_obs sum(r_pdf)]
    fprintf('sum r_pdf is insufficient\n')
    fprintf("imax: %i Rp: %2f p_obs: %2f ks: %2f sumqpdf: %3f\n", i_max, rp, ks, p_obs, sum(r_pdf))
end

%C_sp is probability that a secondary case is classified as a primary
%C_ss is probability that a secondary case is classified as a secondary
if p_obs == 1
    C_sp = 0;
else
    C_sp = sum(r_pdf*p_obs.*exp(log(1-p_obs).*((1:i_max) - 1)));    
end
C_ss = p_obs - C_sp;

%[c,d] = prim_classifier_prob(i_max,rp,rs,ks,p_obs)
%[c,d] = prim_classifier_prob(5,1,0.5,1,0.5)

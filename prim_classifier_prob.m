function [C_pp, C_ps] = prim_classifier_prob(i_max,rp,p_obs)


%p_pdf = probability distribution for primary cases in a cluster
i_arr = 1:i_max;
p_pdf = (1 - 1/rp).^(i_arr-1)/rp;

%q_pdf = probability of a primary case being the ith case in a cluster
temp = 1/rp*cumsum(p_pdf(end:-1:1));
q_pdf = temp(end:-1:1);

if sum(q_pdf) < .99
    % 'sum q_pdf is insufficient'
    % [i_max rp p_obs sum(q_pdf)]
    fprintf('sum q_pdf is insufficient\n')
    fprintf("imax: %i Rp: %2f p_obs: %2f sumqpdf: %3f\n", i_max, rp, p_obs, sum(q_pdf))
end

%C_pp is probability that a primary case is classified as a primary
%C_ps is probability that a primary case is classified as a secondary
if p_obs == 1
    C_pp = q_pdf(1);
else
    C_pp = sum(q_pdf*p_obs.*exp(log(1-p_obs).*((1:i_max) - 1)));    
end
C_ps = p_obs - C_pp;

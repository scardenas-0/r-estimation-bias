function res_arr = sample_classifer( ...
    i_max, ...
    rp_arr, ...
    rs_arr, ...
    ks_arr, ...
    p_obs_arr, ...
    ofile ...
    )

res_arr = zeros(length(rp_arr),length(rs_arr),length(ks_arr),length(p_obs_arr),9);
num_runs = length(res_arr(:))/9;

nn = 0;

% i_max is number of cases to iterate over? highest n
% iterate over arrays for R_p, R_s, k_s, p_obs
for pp = 1:length(rp_arr)
    for ss = 1:length(rs_arr)
        for kk = 1:length(ks_arr)
            for oo = 1:length(p_obs_arr)
                res_arr(pp,ss,kk,oo,:) = ...
                    classifier_performance( ...
                    i_max,rp_arr(pp),rs_arr(ss), ...
                    ks_arr(kk),p_obs_arr(oo));
                nn = nn+1;
                if mod(nn,10) == 0
                    % [nn num_runs]
                    fprintf("%i\t", nn)
                    fprintf("%i\n", num_runs)
                end
            end
        end
    end
end
save(ofile)
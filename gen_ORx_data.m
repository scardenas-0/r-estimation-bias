function gen_ORx_data( ...
    rp_arr, ...
    p_obs_arr, ...
    rs_arr, ...
    inc_arr, ...
    k, ...
    filename ...
    )

res_arr = zeros(length(rp_arr),length(p_obs_arr),length(rs_arr),length(inc_arr),13);
num_runs = length(res_arr(:))/9;

for pp = 1:length(rp_arr)
    for oo = 1:length(p_obs_arr)
        % [rp_arr(pp) p_obs_arr(oo)]
        fprintf("%f\t", rp_arr(pp))
        fprintf("%f\n", p_obs_arr(oo))
        for ss = 1:length(rs_arr)
            for ii = 1:length(inc_arr(:,1))
                res_arr(pp,oo,ss,ii,:) = calc_OR_res(500,rp_arr(pp),rs_arr(ss),k,p_obs_arr(oo),inc_arr(ii,1),inc_arr(ii,2));
            end
        end
    end
end

save(filename)
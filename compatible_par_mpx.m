function  rs_inf_res = compatible_par_mpx(imax,rs_obs_arr,rp,ks,p_obs_arr,ofile)

rs_inf_res = zeros(length(rs_obs_arr),length(p_obs_arr));
for ss = 1:length(rs_obs_arr)
    for pp = 1:length(p_obs_arr)
        [rs_obs_arr(ss) p_obs_arr(pp)]
        rs_mpx_res(ss,pp) = fzero(@(x) ...
            (det_rs_inf(imax, rp,x,ks,p_obs_arr(pp)) ...
            - rs_obs_arr(ss)), 0.5,optimset('TolX',.01));
    end
end

save(ofile)
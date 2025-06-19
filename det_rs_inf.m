function rs_inf = det_rs_inf(i_max,rp,rs,ks,p_obs)

%[rp rs p_obs]
if rs > .80
    rs_inf = 1;
    return
end
if rs < 0
    rs_inf = 0;
    return
end

perf_arr=classifier_performance(i_max,rp,rs,ks,p_obs);
rs_inf = perf_arr(8);

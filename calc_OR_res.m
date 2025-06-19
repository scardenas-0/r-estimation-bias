function OR_res = calc_OR_res(i_max,rp,rs,ks,p_obs,inc_p,inc_s)

c_arr = classifier_performance(i_max,rp,rs,ks,p_obs);
%perf_arr = [C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias];

%true_prim_neg
A = (1-rs)*(1-inc_p);
%true_prim_pos 
B = (1-rs)*(inc_p);
%true_sec_neg 
C= (rs)*(1-inc_s);
%true_sec_pos 
D = (rs)*(inc_s);

%obs_prim_neg
Ax = A*c_arr(1) + C*c_arr(3);
%obs_prim_pos
Bx = B*c_arr(1) + D*c_arr(3);
%obs_sec_neg
Cx = A*c_arr(2) + C*c_arr(4);
%obs_sec_pos
Dx = B*c_arr(2) + D*c_arr(4);

% OR = inc_p*(1-inc_s)/((1-inc_p)*inc_s)
OR = B*C/(A*D);
ORx = Bx*Cx/(Ax*Dx);

%Observed primary incidence
inc_px = Bx/(Ax+Bx);
inc_sx = Dx/(Cx+Dx);

OR_res = [OR ORx inc_px inc_sx c_arr];

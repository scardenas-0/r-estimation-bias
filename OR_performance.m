function res_arr = OR_performance( ...
    p_obs, ...
    C_pp, C_ps, C_sp, C_ss, ...
    X_pn, X_pp, X_sn, X_sp ...
    )
[C_pp C_ps C_sp C_ss]
% 84/248, 85/248, 74/248, 5/248

%[p_obs C_pp C_ps C_sp C_ss X_pn X_pp X_sn X_sp]
Q_pn = p_obs * (X_sn*C_sp - X_pn*C_ss)/(C_ps*C_sp - C_pp*C_ss);
Q_pp = p_obs * (X_sp*C_sp - X_pp*C_ss)/(C_ps*C_sp - C_pp*C_ss);
Q_sn = p_obs * (X_sn - X_pn*C_ps/C_pp)/(C_ss - C_sp*C_ps/C_pp);
Q_sp = p_obs * (X_sp - X_pp*C_ps/C_pp)/(C_ss - C_sp*C_ps/C_pp); % this one is negative for some reason
% Q_sp_int_1 = X_sp;
% Q_sp_int_2 = X_pp*C_ps/C_pp;
% Q_sp_int_1
% Q_sp_int_2
Q_sp
res_arr = [Q_pn Q_pp Q_sn Q_sp];

i_max = 800;
ks_arr = [0.3 1 1e3];
p_obs_arr = [0.1:0.05:1];

%Theory
res_arr = sample_classifer( ...
    i_max, ...
    1:.05:2, ...
    0.2:0.3:0.8, ...
    ks_arr, ...
    p_obs_arr, ...
    'Data/011119_rp_pobs' ...
    );

%%
% Observed OR example
clear;
k = 0.3;
rp_arr = 1:.05:2;
p_obs_arr = .25:.05:1;
rs_arr = [0.25, 0.5 0.75];
inc_arr = [0.8 0.5; 0.5 0.2; 0.1 0.027027];
gen_ORx_data( ...
    rp_arr, ...
    p_obs_arr, ...
    rs_arr, ...
    inc_arr, ...
    k, ...
    'Data/040419_ORx_example' ...
    )

%%
% Measles
% https://www.cdc.gov/measles/downloads/report-elimination-measles-rubella-crs.pdf
% Data from Table 5 (with 'all chains of transmission' noted to be
% incorrectly summed in report.
% Total cases 904 cases among 444 chains
% Compat_arr = compatible_par_inf(400,1-444/904,1:.05:1.5,1,0.5:.05:1,'Data/010919_msls_model');
% New measles data 
% in 2020-2024: cases = 338, chains = 92
% in 2020-2022: cases = 183, chains = 31
% in 2023-2024: cases = 155, chains = 61
% in 2020-2021: cases = 62, chains = 16
% in 2022-2024: cases = 276, chains = 76
clear;
k_msls = 0.3;

% Compat_arr = compatible_par_inf(800,1-444/904,1:.25:1.3,k_msls,0.25:.25:1,'Data/040419_msls_model');
% Compat_arr = compatible_par_inf(800,1-444/904,1:.025:1.3,k_msls,0.25:.025:1,'Data/040419_msls_model');
% Compat_arr_2024 = compatible_par_inf(800, 1-92/338, 1:.05:1.5, k_msls, 0.25:.05:1, 'Data/031125_msls_model');
Compat_arr_2024 = compatible_par_inf( ...
    800, ...
    1-31/183, ...
    1:.05:1.5, ...
    k_msls, ...
    0.25:.05:1, ...
    'Data/061925_msls_02_model' ...
    );

% d_msls = load('Data/040419_msls_r');
d_msls = load('Data/061925_msls_02_model');
rp_num = length(d_msls.rp_arr);
po_num = length(d_msls.p_obs_arr);
res_arr = zeros(rp_num,po_num,9);

for pp = 1:rp_num
    for oo = 1:po_num
        [pp rp_num oo po_num]
        res_arr(pp,oo,:) = classifier_performance( ...
            800, ...
            d_msls.rp_arr(pp), ...
            d_msls.rs_inf_res(pp,oo), ...
            k_msls, ...
            d_msls.p_obs_arr(oo) ...
            );
    end
end
save('Data/061925_msls_02_perf')

%%

% Monkeypox
% PNAS article 760 cases
% Ecohealth article 156 sites
% Compat_arr = compatible_par_inf(500,604/760,1:.1:2,0.3,0.9:.01:1,'Data/012518_mpx_model');
%
% Re-framing it so that we might observe between 200 and 700 infection
% clusters
%rs_inf_res = compatible_par_inf(500,604/760,1:.2:4,0.3,0.25:.025:1,'Data/011219_mpx_model');
%compat_arr = compatible_par_mpx(500,0.2:.025:0.8,1,0.3,0.25:.025:1,'Data/011119_mpx_model');
%rs_inf_res = compatible_par_inf(imax,rs_obs,rp_arr,ks,p_obs_arr,ofile)

% Approximate numbers for MRSA
% Based on Shea2018_phillips_final.pptx, slide #25 (OLD)
% There were 392 sequences 42 'transmission events' involving 96 patients.
% Will interpret that as meaning there were 392-96+42 = 338 primary
% infections and 96-42 = 54 secondary infections
% OLD - not usinganymore: compatible_par_inf(500,54/392,.8:.025:1.1,0.3,0.5:.05:1,'Data/050918_mrsa_model');
%
% MRSA script
%MRSA_example_122018

%%
% MPX 1980s
% Primary infection: Male 142 Female 103
% Secondary infection: Male 40 Female 53

%%
% Nipah virus 2001-2014
% Primary infection: Male 84 Female 85
% Secondary infection: Male 74 Female 5

clear;
% k_mpx = 0.3
k_nipah = 0.3
% Compat_arr = compatible_par_inf(800,93/338,1:.025:1.3,k_mpx,0.25:.025:1,'Data/040419_mpx_r');
Compat_arr = compatible_par_inf( ...
    800, ...
    79/248, ...
    1:.005:1.1, ...
    k_nipah, ...
    0.25:.025:1, ...
    'Data/050125_nipah_r' ...
    );
%%
% MPX 1980s part 2

clear
% d_mpx = load('Data/040419_mpx_r')
d_nipah = load('Data/050125_nipah_r')

rp_num = length(d_nipah.rp_arr);
po_num = length(d_nipah.p_obs_arr);
res_arr = zeros(rp_num,po_num,13);
k_nipah = 0.3
for pp = 1:rp_num
    for oo = 1:po_num
        [pp rp_num oo po_num]
        res_arr(pp,oo,1:9) = classifier_performance( ...
            800, d_nipah.rp_arr(pp), d_nipah.rs_inf_res(pp,oo), k_nipah, d_nipah.p_obs_arr(oo) ...
            );
        % Treating female as 'positive'
        res_arr(pp,oo,10:13) = OR_performance( ...
            d_nipah.p_obs_arr(oo), ...
            res_arr(pp,oo,1), res_arr(pp,oo,2), res_arr(pp,oo,3), res_arr(pp,oo,4), ...
            84/248, 85/248, 74/248, 5/248 ...
            );
            % measles: 142/338, 103/338, 40/338, 53/338 ...
            % nipah: 84/248, 85/248, 74/248, 5/248
        % res_arr(pp,oo,13)
    end
end

save('Data/050125_nipah_OR')

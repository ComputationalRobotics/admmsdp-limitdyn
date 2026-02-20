clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

addpath(genpath("../utils/"));
addpath(genpath("../admm-single-light/"));

%% input and output
fasrc_path = "/n/netscratch/hankyang_lab/Lab/ksc19/2025-field/codes/";
if exist(fasrc_path, 'dir')
    input_path = fasrc_path + "data/";
    output_path = fasrc_path + "sig-wave/data/";
else
    input_path = "../data/";
    output_path = "./data/";
end
type = "neu2g";
SDP_data = load(input_path + type + "/SDP_data.mat");
if isfield(SDP_data, 'SDP_data')
    SDP_data = SDP_data.SDP_data;
end
SDP_data.A = SDP_data.At';
if ~exist(output_path + type, "dir")
    mkdir(output_path + type);
end
% set up logs
if ~exist("./logs/" + type, "dir")
    mkdir("./logs/" + type);
end
diary("./logs/" + type + "/generate_start_points.txt");
diary on;

fprintf("generating start points... type: %s \n\n", type);

%% gpu-specific
[gpu_ok, ~] = if_has_gpu();
vec_len = size(SDP_data.C, 1);
mat_size = (-1 + sqrt(1 + 8 * vec_len)) / 2;
% if mat_size > 5,000, consider gpu for projection
if gpu_ok && mat_size > 5000
    SDP_data.if_gpu = true;
    fprintf("use GPU for PSD projection! \n");
end

%% all-zero initial guess
initial_guess_info.X_mat = zeros(size(smat_single(SDP_data.C), 1)); 
initial_guess_info.X_mat = 0.5 * (initial_guess_info.X_mat + initial_guess_info.X_mat');
initial_guess_info.S_mat = zeros(size(initial_guess_info.X_mat));
initial_guess_info.S_mat = 0.5 * (initial_guess_info.S_mat + initial_guess_info.S_mat');
sig_data.sig = 1.0;
Z_init = svec_single(initial_guess_info.X_mat - sig_data.sig * initial_guess_info.S_mat);

%% using update_sig_1
sig_data = struct();
sig_data.stop_iter = 2e4;
sig_data = set_sig(1.0, sig_data);
% initialize
iter_data = set_iter(Z_init, sig_data);
[iter_data_prev, SDP_data] = scale_data(iter_data, sig_data, SDP_data, true); 
iter_data = three_step_ADMM(iter_data_prev, sig_data, SDP_data); 
iter_data_next = three_step_ADMM(iter_data, sig_data, SDP_data); 
log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
draw_data = set_draw(0, [], sig_data, log_data, true);
% start, after stop_iter, run another stop_iter iteration for stabilization
tol = 1e-10;
for iter = 1: sig_data.stop_iter * 2
    iter_data_prev = iter_data;
    iter_data = iter_data_next;
    iter_data_next = three_step_ADMM(iter_data_next, sig_data, SDP_data);
    log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
    draw_data = set_draw(iter, draw_data, sig_data, log_data);
    sig_data = update_sig_1(iter_data.cnt, sig_data, log_data);

    if mod(iter, 10) == 1 
        fprintf("iter: %d, [%3.2e, %3.2e, %d], [%3.2e | %3.2e | %3.2e], [%3.2e, %3.2e], %3.2e \n",...
                    iter, log_data.Z_diff_norm, log_data.Z_diff_ang, log_data.Z_rank, ...
                    log_data.pinf, log_data.dinf, log_data.relgap, ...
                    log_data.pobj, log_data.dobj, sig_data.sig);
    end

    if log_data.kkt < tol
        fprintf("Convergent! \n");
        break;
    end
end
% save 
iter_data_scaleback = scaleback_data(iter_data, sig_data, SDP_data);
iter_data_scaleback.sig_data = sig_data;
save(output_path + type + "/iter_data_start.mat", "iter_data_scaleback");
save(output_path + type + "/draw_data_start.mat", "draw_data");

diary off;




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
type = "swissroll";
SDP_data = load(input_path + type + "/SDP_data.mat");
if isfield(SDP_data, 'SDP_data')
    SDP_data = SDP_data.SDP_data;
end
SDP_data.A = SDP_data.At';
output_data_path = output_path + type + "/sig-wave/";
output_figs_path = "./figs/" + type + "/sig-wave/";
if ~exist(output_data_path, "dir")
    mkdir(output_data_path);
else
    empty_folder(output_data_path);
end
% set up logs
if ~exist("./logs/" + type, "dir")
    mkdir("./logs/" + type);
end
diary("./logs/" + type + "/study_sig_wave.txt");
diary on;

fprintf("study sig wave... type: %s \n\n", type);

%% gpu-specific
[gpu_ok, ~] = if_has_gpu();
vec_len = size(SDP_data.C, 1);
mat_size = (-1 + sqrt(1 + 8 * vec_len)) / 2;
% if mat_size > 5,000, consider gpu for projection
if gpu_ok && mat_size > 5000
    SDP_data.if_gpu = true;
    fprintf("use GPU for PSD projection! \n");
end

%% load initial iter_data
load(output_path + type + "/iter_data_start.mat");
sig_data = iter_data_scaleback.sig_data;
[iter_data_base, SDP_data] = scale_data(iter_data_scaleback, sig_data, SDP_data, true);
iter_data = three_step_ADMM(iter_data_base, sig_data, SDP_data);
% for anchor angle test
X_diff_anchor = iter_data.X - iter_data_base.X;
S_diff_anchor = iter_data.S - iter_data_base.S;

%% test sig wave
sig_data.a = 0;
sig_data.b = 1/4;
sig_data.sig_base = sig_data.sig;
sig_data.N = ceil(2e4);
Ampl_list = [2, 1e1, 1e2];

for ii = 1: length(Ampl_list)
    fprintf("\n --- Ampl = %3.2e, N = %d --- \n", Ampl_list(ii), sig_data.N);
    sig_data.Ampl = Ampl_list(ii);

    % forward pass
    sig_data.sig = sig_data.sig_base;
    iter_data_prev = iter_data_base;
    iter_data = three_step_ADMM(iter_data_prev, sig_data, SDP_data); 
    iter_data_next = three_step_ADMM(iter_data, sig_data, SDP_data); 
    log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
    draw_data_forward = set_draw(0, [], sig_data, log_data, true);
    X_anchor_diff_ang_list = [];
    S_anchor_diff_ang_list = [];
    for iter = 1: ceil(sig_data.N/4)
        iter_data_prev = iter_data;
        iter_data = iter_data_next;
        iter_data_next = three_step_ADMM(iter_data_next, sig_data, SDP_data);
        log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
        draw_data_forward = set_draw(iter, draw_data_forward, sig_data, log_data);
        sig_data = update_sig_wave(iter, sig_data, log_data);
        
        X_diff = iter_data_next.X - iter_data.X;
        S_diff = iter_data_next.S - iter_data.S;
        X_anchor_diff_ang_list = [X_anchor_diff_ang_list, get_ang(X_diff_anchor, X_diff)];
        S_anchor_diff_ang_list = [S_anchor_diff_ang_list, get_ang(S_diff_anchor, S_diff)];
    
        if mod(iter, 10) == 1 
            fprintf("iter: %d, [%3.2e, %3.2e, %d], [%3.2e | %3.2e | %3.2e], [%3.2e, %3.2e], %3.2e \n",...
                        iter, log_data.Z_diff_norm, log_data.Z_diff_ang, log_data.Z_rank, ...
                        log_data.pinf, log_data.dinf, log_data.relgap, ...
                        log_data.pobj, log_data.dobj, sig_data.sig);
        end
    end
    draw_data_forward.X_anchor_diff_ang_list = X_anchor_diff_ang_list;
    draw_data_forward.S_anchor_diff_ang_list = S_anchor_diff_ang_list;

    % backward pass
    sig_data.sig = sig_data.sig_base;
    iter_data_prev = iter_data_base;
    iter_data = three_step_ADMM(iter_data_prev, sig_data, SDP_data); 
    iter_data_next = three_step_ADMM(iter_data, sig_data, SDP_data); 
    log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
    draw_data_backward = set_draw(0, [], sig_data, log_data, true);
    X_anchor_diff_ang_list = [];
    S_anchor_diff_ang_list = [];
    for iter = ceil(2*sig_data.N/4): ceil(3*sig_data.N/4)
        iter_data_prev = iter_data;
        iter_data = iter_data_next;
        iter_data_next = three_step_ADMM(iter_data_next, sig_data, SDP_data);
        log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
        draw_data_backward = set_draw(iter, draw_data_backward, sig_data, log_data);
        sig_data = update_sig_wave(iter, sig_data, log_data);
        
        X_diff = iter_data_next.X - iter_data.X;
        S_diff = iter_data_next.S - iter_data.S;
        X_anchor_diff_ang_list = [X_anchor_diff_ang_list, get_ang(X_diff_anchor, X_diff)];
        S_anchor_diff_ang_list = [S_anchor_diff_ang_list, get_ang(S_diff_anchor, S_diff)];
    
        if mod(iter, 10) == 1 
            fprintf("iter: %d, [%3.2e, %3.2e, %d], [%3.2e | %3.2e | %3.2e], [%3.2e, %3.2e], %3.2e \n",...
                        iter, log_data.Z_diff_norm, log_data.Z_diff_ang, log_data.Z_rank, ...
                        log_data.pinf, log_data.dinf, log_data.relgap, ...
                        log_data.pobj, log_data.dobj, sig_data.sig);
        end
    end
    draw_data_backward.X_anchor_diff_ang_list = X_anchor_diff_ang_list;
    draw_data_backward.S_anchor_diff_ang_list = S_anchor_diff_ang_list;

    % save data
    draw_data.forward = draw_data_forward;
    draw_data.backward = draw_data_backward;
    draw_data_filename = output_data_path + sprintf("draw_data_Ampl=%d_N=%d.mat", ceil(sig_data.Ampl), sig_data.N);
    save(draw_data_filename, "draw_data");
end

diary off;

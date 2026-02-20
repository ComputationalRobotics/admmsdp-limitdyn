clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

addpath(genpath("../utils/"));
addpath(genpath("../admm-single-light/"));

%% hyper-parameters
type = "swissroll";
t = datetime('now','Format','yyyyMMdd-HHmmss');  
curr_time = char(t);
% for warm start usage only
% if warmstart, only support "baseline", "linesearch", and "sigwave"
warmstart_type = "";
last_time = "";
last_iter = 0;
if ~any(warmstart_type == ["", "baseline", "linesearch", "sigwave"])
    error("unsupported warmstart type! \n");
end
% solving related 
maxiter = ceil(1e6);
tol = 1e-10;
save_interval = ceil(1e5);

%% set up logs
log_path = "./logs/" + type + "/" + "baseline/" + curr_time;
if ~exist(log_path, "dir")
    mkdir(log_path);
end
diary(log_path + "/main_baseline.txt");
diary on;

%% input and output
fasrc_path = "/n/netscratch/hankyang_lab/Lab/ksc19/2025-field/codes/";
if exist(fasrc_path, 'dir')
    input_path = fasrc_path + "acceleration/data/" + type + "/";
    input_SDP_path = fasrc_path + "data/" + type + "/";
    input_iter_path = fasrc_path + "sig-wave/data/" + type + "/";
    output_iter_path = fasrc_path + "acceleration/data/" + type + "/";
else
    input_path = "./data/" + type + "/";
    input_SDP_path = "../data/" + type + "/";
    input_iter_path = "../sig-wave/data/" + type + "/";
    output_iter_path = "./data/" + type + "/";
end
SDP_data = load(input_SDP_path + "/SDP_data.mat");
if isfield(SDP_data, 'SDP_data')
    SDP_data = SDP_data.SDP_data;
end
SDP_data.A = SDP_data.At';
% update input_iter_path and load input_iter
if warmstart_type ~= "" && last_time ~= ""
    input_iter_path = input_path + warmstart_type + "/" + last_time + "/";
    filename = input_iter_path + "iter_data_" + string(last_iter) + ".mat";
    if ~isfile(filename)
        error("no such warmstart file! \n");
    else 
        fprintf("warmstart! loading %s... \n", filename);
        load(filename);
    end
else 
    filename = input_iter_path + "iter_data_start.mat";
    fprintf("no warmstart! loading %s... \n", filename);
    load(filename);
end
% update output_iter_path
output_iter_path = output_iter_path + "/baseline/" + curr_time + "/";
if ~exist(output_iter_path, "dir")
    mkdir(output_iter_path);
end

fprintf("main baseline... type: %s \n\n", type);

%% gpu-specific
[gpu_ok, ~] = if_has_gpu();
vec_len = size(SDP_data.C, 1);
mat_size = (-1 + sqrt(1 + 8 * vec_len)) / 2;
% if mat_size > 5,000, consider gpu for projection
if gpu_ok && mat_size > 5000
    SDP_data.if_gpu = true;
    fprintf("use GPU for PSD projection! \n");
end

tic;

%% SDP solving
% initialize
sig_data = iter_data_scaleback.sig_data;
[iter_data_prev, SDP_data] = scale_data(iter_data_scaleback, sig_data, SDP_data, true); 
iter_data = three_step_ADMM(iter_data_prev, sig_data, SDP_data); 
iter_data_next = three_step_ADMM(iter_data, sig_data, SDP_data); 
log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
draw_data = set_draw(0, [], sig_data, log_data, true);
% solve
for iter = 1: maxiter
    iter_data_prev = iter_data;
    iter_data = iter_data_next;
    iter_data_next = three_step_ADMM(iter_data_next, sig_data, SDP_data);
    log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
    draw_data = set_draw(iter, draw_data, sig_data, log_data);

    if mod(iter, 50) == 1 
        fprintf("iter: %d, [%3.2e, %3.2e, %d], [%3.2e | %3.2e | %3.2e], [%3.2e, %3.2e], %3.2e \n",...
                    iter, log_data.Z_diff_norm, log_data.Z_diff_ang, log_data.Z_rank, ...
                    log_data.pinf, log_data.dinf, log_data.relgap, ...
                    log_data.pobj, log_data.dobj, sig_data.sig);
    end

    if log_data.kkt < tol
        fprintf("Convergent! \n");
        break;
    end

    % regularly save files
    if mod(iter, save_interval) == 0
        fprintf("saving iter_data and draw_data in iter: %d... \n", iter);
        iter_data_scaleback = scaleback_data(iter_data, sig_data, SDP_data);
        iter_data_scaleback.sig_data = sig_data;
        iter_data_scaleback.warmstart_type = warmstart_type;
        iter_data_scaleback.last_time = last_time;
        iter_data_scaleback.last_iter = last_iter;
        save(output_iter_path + "iter_data_" + string(iter) + ".mat", "iter_data_scaleback");
        save(output_iter_path + "draw_data_" + string(iter) + ".mat", "draw_data");
    end
end

elapsedTime = toc;   
fprintf('Elapsed time: %.4f seconds\n', elapsedTime);

%% save final data
fprintf("saving iter_data and draw_data in iter: %d... \n", iter);
iter_data_scaleback = scaleback_data(iter_data, sig_data, SDP_data);
iter_data_scaleback.sig_data = sig_data;
iter_data_scaleback.warmstart_type = warmstart_type;
iter_data_scaleback.last_time = last_time;
iter_data_scaleback.last_iter = last_iter;
save(output_iter_path + "iter_data_" + string(iter) + ".mat", "iter_data_scaleback");
save(output_iter_path + "draw_data_" + string(iter) + ".mat", "draw_data");

diary off;
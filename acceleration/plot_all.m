clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

addpath(genpath("../utils/"));
addpath(genpath("../admm-single-light/"));

%% type
% supported type: 1dc1024, 1et2048, 1tc2048, 1zc1024, cancer100, cnhil10, cphil12, foot, G40mb, G48mb, G48mc, 
% hamming8, hamming9, hand, neosfbr25, neosfbr30e8, neu1g, neu2g, neu3g, r12000, swissroll, theta12, theta102,
% theta123, texture
type_list = [
    "1dc1024", "1et2048", "1tc2048", "1zc1024", "cancer100", "cnhil10", "cphil12", "foot", "G40mb", "G48mb", "G48mc", "hamming8", "hamming9", "hand", "neosfbr25", "neosfbr30e8", "neu1g", "neu2g", "neu3g", "r12000", "swissroll", "theta12", "theta102", "theta123", "texture"
];

for type = type_list

    %% input and output
    fasrc_path = "/n/netscratch/hankyang_lab/Lab/ksc19/2025-field/codes/";
    if exist(fasrc_path, 'dir')
        input_data_path = fasrc_path + "acceleration/data/" + type + "/";
        input_startpoint_path = fasrc_path + "sig-wave/data/" + type + "/";
    else
        input_data_path = "./data/" + type + "/";
        input_startpoint_path = "../sig-wave/data/" + type + "/";
    end
    % figs output
    output_figs_path = "./figs/" + type + "/";
    if ~exist(output_figs_path, 'dir') 
        mkdir(output_figs_path);
    end

    %% load data
    if_bas = false;
    if_lin = false;
    if_sig = false;
    % first load startpoint data
    draw_data_sta = load(input_startpoint_path + "draw_data_start.mat");
    draw_data_sta = draw_data_sta.draw_data;
    % load baseline, linesearch, sigwave data
    filenames = get_filename_type(type);
    if ~isempty(filenames.baseline)
        if_bas = true;
        draw_data_bas = load(input_data_path + filenames.baseline.draw_data_filename); 
        draw_data_bas = draw_data_bas.draw_data;
    else
        draw_data_bas = [];
    end
    if ~isempty(filenames.linesearch)
        if_lin = true;
        draw_data_lin = load(input_data_path + filenames.linesearch.draw_data_filename); 
        draw_data_lin = draw_data_lin.draw_data;
    else 
        draw_data_lin = [];
    end
    if ~isempty(filenames.sigwave)
        if_sig = true;
        draw_data_sig = load(input_data_path + filenames.sigwave.draw_data_filename); 
        draw_data_sig = draw_data_sig.draw_data;
    else 
        draw_data_sig = [];
    end

    delete(output_figs_path + "log.txt");
    diary(output_figs_path + "log.txt");
    diary on;

    %% plot normal kkt, normal Z_diff_ang, normal Z_diff_norm
    draw_original(draw_data_sta, draw_data_bas, output_figs_path);
    continue;

    %% plot linesearch and baseline's kkt and Znorm
    if ~isempty(draw_data_lin)
        draw_linesearch(draw_data_bas, draw_data_lin, output_figs_path);
    end

    %% plot sigwave and baseline's kkt, Xdiffnorm, Sdiffnorm, Xdiffang, Sdiffang
    if ~isempty(draw_data_sig)
        draw_sigwave(draw_data_bas, draw_data_sig, output_figs_path);
    end

    diary off;

end


%% helper functions
function filenames = get_filename_type(type)
    filenames = struct();
    switch type
        case "1dc1024"
            filenames.baseline   = get_info("baseline",   "20251222-151734", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-201352", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-193342", 1000000);
        case "1et2048"
            filenames.baseline   = get_info("baseline",   "20251225-140541", 69237);
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "1tc2048"
            filenames.baseline   = get_info("baseline",   "20251225-144610", 1000000);
            filenames.linesearch = get_info("linesearch", "20251226-132856", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251226-133024", 1000000);
        case "1zc1024"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "cancer100"
            filenames.baseline   = get_info("baseline",   "20251222-154206", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-195442", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-192017", 1000000);
        case "cnhil10"
            filenames.baseline   = get_info("baseline",   "20251222-144657", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-195217", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-154208", 1000000);
        case "cphil12"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "foot"
            filenames.baseline   = get_info("baseline",   "20251222-154938", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-202622", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-194559", 1000000);
        case "G40mb"
            filenames.baseline   = get_info("baseline",   "20251225-155604", 1000000);
            filenames.linesearch = get_info("linesearch", "20251226-133204", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251226-133318", 1000000);
        case "G48mb"
            filenames.baseline   = get_info("baseline",   "20251222-152646", 72379);
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "G48mc"
            filenames.baseline   = get_info("baseline",   "20251225-160059", 91203);
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "hamming8"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "hamming9"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "hand"
            filenames.baseline   = get_info("baseline",   "20251222-160203", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-203307", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-195330", 1000000);
        case "neosfbr25"
            filenames.baseline   = get_info("baseline",   "20251222-161110", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-204010", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-200723", 1000000);
        case "neosfbr30e8"
            filenames.baseline   = get_info("baseline",   "20251222-162000", 800000);
            filenames.linesearch = get_info("linesearch", "20251223-204600", 800000);
            filenames.sigwave    = get_info("sigwave",    "20251224-201116", 700000);
        case "neu1g"
            filenames.baseline   = get_info("baseline",   "20251222-162031", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-205156", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-185829", 80378);
        case "neu2g"
            filenames.baseline   = get_info("baseline",   "20251222-130424", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-195156", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-154408", 1000000);
        case "neu3g"
            filenames.baseline   = get_info("baseline",   "20251222-165918", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-205457", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-190427", 1000000);
        case "r12000"
            filenames.baseline   = get_info("baseline",   "20251225-160159", 900000);
            filenames.linesearch = get_info("linesearch", "20251226-143018", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251226-173447", 1000000);
        case "swissroll"
            filenames.baseline   = get_info("baseline",   "20251225-160229", 1000000);
            filenames.linesearch = get_info("linesearch", "20251226-191025", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251226-193633", 1000000);
        case "theta12"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "theta102"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "theta123"
            filenames.baseline   = [];
            filenames.linesearch = [];
            filenames.sigwave    = [];
        case "texture"
            filenames.baseline   = get_info("baseline",   "20251222-162723", 1000000);
            filenames.linesearch = get_info("linesearch", "20251223-210307", 1000000);
            filenames.sigwave    = get_info("sigwave",    "20251224-201950", 1000000);
        otherwise
            error("dataset type not supported! \n");
    end
end

function info = get_info(technique, date, iter)
    info.technique = technique;
    info.date = date;
    info.iter = iter;
    info.draw_data_filename = technique + "/" + date + "/" + "draw_data_" + string(iter) + ".mat";
    info.iter_data_filename = technique + "/" + date + "/" + "iter_data_" + string(iter) + ".mat";
end

%% functions for draw original
function draw_original(draw_data_sta, draw_data_bas, img_prefix)
    % Data
    if isempty(draw_data_bas)
        Z_diff_norm_list = draw_data_sta.Z_diff_norm_list;
        Z_diff_ang_list  = draw_data_sta.Z_diff_ang_list;
        pinf_list        = draw_data_sta.pinf_list;
        dinf_list        = draw_data_sta.dinf_list;
        relgap_list      = draw_data_sta.relgap_list;
        bas_pobj_list    = [];
        bas_dobj_list    = [];
    else 
        Z_diff_norm_list = [draw_data_sta.Z_diff_norm_list, draw_data_bas.Z_diff_norm_list];
        Z_diff_ang_list  = [draw_data_sta.Z_diff_ang_list, draw_data_bas.Z_diff_ang_list];
        pinf_list        = [draw_data_sta.pinf_list, draw_data_bas.pinf_list];
        dinf_list        = [draw_data_sta.dinf_list, draw_data_bas.dinf_list];
        relgap_list      = [draw_data_sta.relgap_list, draw_data_bas.relgap_list];
        bas_pobj_list    = draw_data_bas.pobj_list;
        bas_dobj_list    = draw_data_bas.dobj_list;
    end
    kkt_list = max([pinf_list; dinf_list; relgap_list]);

    % first plot I: original data's kkt 
    rows = 150; cols = 300; linewidth = 1.2; ms = 12; xls = 12; yls = 12; fontsize = 12;
    LineColors = linspecer(3, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % semilogy plot, two lines with different colors
    plot(log10(kkt_list), 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on;
    plot(log10(Z_diff_norm_list), 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    grid on;
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f1, img_prefix + "original_kkt.png", '-dpng', ['-r' num2str(500)]); close(f1);

    % first plot II: original data's ang
    rows = 150; cols = 300; linewidth = 1.2; ms = 12; xls = 12; yls = 12; fontsize = 12;
    LineColors = linspecer(3, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % semilogy plot, two lines with different colors
    plot(log10(Z_diff_ang_list), 'LineWidth', linewidth, 'Color', LineColors(3, :)); grid on;
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f1, img_prefix + "original_ang.png", '-dpng', ['-r' num2str(500)]); close(f1);

    % second plot: original legend
    rows = 50; cols = 600; linewidth = 3;
    LineColors = linspecer(3, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f2 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    legend_1 = "$\log_{10}(r_{\max}^{(k)})$"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$\log_{10}(\parallel \Delta Z^{(k)} \parallel_F)$"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend_3 = "$\log_{10}(\angle (\Delta Z^{(k)}, \Delta Z^{(k+1)}))$"; plot(nan, 'Color', LineColors(3, :), 'LineWidth', linewidth); 
    legend({legend_2, legend_1, legend_3}, 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f2, img_prefix + "original_legend.png", '-dpng', ['-r' num2str(500)]); close(f2);

    if isempty(draw_data_bas)
        return;
    end

    % return directly
    return;

    % third plot: pobj and dobj data
    rows = 300; cols = 300; linewidth = 1.2; 
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f3 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % semilogy plot, two lines with different colors
    semilogy(bas_pobj_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    semilogy(bas_dobj_list, 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on;
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex'); grid on;
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    % Save & close
    print(f3, img_prefix + "baseline_obj_data.png", '-dpng', ['-r' num2str(500)]); close(f3);

    % fourth plot: original legend
    rows = 50; cols = 500; linewidth = 3;
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f4 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    legend_1 = "$p^{(k)}$"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$d^{(k)}$"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend({legend_1, legend_2}, 'Interpreter', 'latex', 'NumColumns', 2, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f4, img_prefix + "baseline_obj_legend.png", '-dpng', ['-r' num2str(500)]); close(f4);
end

%% functions for linesearch
function draw_linesearch(draw_data_bas, draw_data_lin, img_prefix)
    % Data
    bas_kkt_list = max([draw_data_bas.pinf_list; draw_data_bas.dinf_list; draw_data_bas.relgap_list]);
    lin_kkt_list = max([draw_data_lin.pinf_list; draw_data_lin.dinf_list; draw_data_lin.relgap_list]);
    bas_Z_norm_list = draw_data_bas.Z_norm_list;
    lin_Z_norm_list = draw_data_lin.Z_norm_list;
    len = min([length(bas_kkt_list), length(lin_kkt_list)]) - 1;
    bas_kkt_list = bas_kkt_list(1: len); lin_kkt_list = lin_kkt_list(1: len);
    bas_Z_norm_list = bas_Z_norm_list(1: len); lin_Z_norm_list = lin_Z_norm_list(1: len);

    % print info
    [best_kkt, best_id] = min(bas_kkt_list);
    fprintf("baseline best kkt: %3.2e, reached at iteration %d \n", best_kkt, best_id);
    [best_kkt, best_id] = min(lin_kkt_list);
    fprintf("linesearch best kkt: %3.2e, reached at iteration %d \n", best_kkt, best_id);
    if any(lin_kkt_list(end-100:end) > 1e-1)
        fprintf("linesearch diverging... \n");
    end

    % first plot: kkt compare data
    rows = 300; cols = 300; linewidth = 1.2; ms = 12; xls = 13; yls = 13; fontsize = 13;
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % semilogy plot, two lines with different colors
    semilogy(bas_kkt_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on; 
    semilogy(lin_kkt_list, 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on; grid on;
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    % Save & close
    print(f1, img_prefix + "linesearch_kkt_comp_data.png", '-dpng', ['-r' num2str(500)]); close(f1);

    % second plot: Znorm compare data
    % Prepare figure
    f2 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % semilogy plot, two lines with different colors
    semilogy(bas_Z_norm_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on; 
    semilogy(lin_Z_norm_list, 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on; grid on;
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    % Save & close
    print(f2, img_prefix + "linesearch_Znorm_comp_data.png", '-dpng', ['-r' num2str(500)]); close(f2);

    % third plot: kkt compare legend
    rows = 50; cols = 500; linewidth = 3;
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f3 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    legend_1 = "$r_{\max}$: Baseline"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$r_{\max}$: Line Search"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend({legend_1, legend_2}, 'Interpreter', 'latex', 'NumColumns', 2, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f3, img_prefix + "linesearch_kkt_comp_legend.png", '-dpng', ['-r' num2str(500)]); close(f3);

    % fourth plot: Znorm compare legend
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f4 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    legend_1 = "$\parallel Z^{(k)} \parallel_F$: Baseline"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$\parallel Z^{(k)} \parallel_F$: Line Search"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend({legend_1, legend_2}, 'Interpreter', 'latex', 'NumColumns', 2, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f4, img_prefix + "linesearch_Znorm_comp_legend.png", '-dpng', ['-r' num2str(500)]); close(f4);
end

%% functions for sigwave
function draw_sigwave(draw_data_bas, draw_data_sig, img_prefix)
    % Data
    bas_kkt_list         = max([draw_data_bas.pinf_list; draw_data_bas.dinf_list; draw_data_bas.relgap_list]);
    sig_kkt_list         = max([draw_data_sig.pinf_list; draw_data_sig.dinf_list; draw_data_sig.relgap_list]);
    sig_X_diff_norm_list = draw_data_sig.X_diff_norm_list;
    sig_S_diff_norm_list = draw_data_sig.S_diff_norm_list;
    sig_X_diff_ang_list  = draw_data_sig.X_diff_ang_list;
    sig_S_diff_ang_list  = draw_data_sig.S_diff_ang_list;
    sig_sig_list         = draw_data_sig.sig_list;
    len = min([length(bas_kkt_list), length(sig_kkt_list)]) - 1;
    bas_kkt_list = bas_kkt_list(1: len); sig_kkt_list = sig_kkt_list(1: len);
    sig_X_diff_norm_list = sig_X_diff_norm_list(1: len); sig_S_diff_norm_list = sig_S_diff_norm_list(1: len);
    sig_X_diff_ang_list = sig_X_diff_ang_list(1: len); sig_S_diff_ang_list = sig_S_diff_ang_list(1: len);
    sig_sig_list = sig_sig_list(1: len);

    % print info
    [best_kkt, best_id] = min(bas_kkt_list);
    fprintf("baseline best kkt: %3.2e, reached at iteration %d \n", best_kkt, best_id);
    [best_kkt, best_id] = min(sig_kkt_list);
    fprintf("sigwave best kkt: %3.2e, reached at iteration %d \n", best_kkt, best_id);
    if any(sig_kkt_list(end-100:end) > 1e-1)
        fprintf("sigwave diverging... \n");
    end

    % first plot: kkt compare data
    rows = 300; cols = 300; linewidth = 1.2; ms = 12; xls = 13; yls = 13; fontsize = 13;
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % semilogy plot, two lines with different colors
    semilogy(bas_kkt_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on; 
    semilogy(sig_kkt_list, 'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on; grid on;
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    % Save & close
    print(f1, img_prefix + "sigwave_kkt_comp_data.png", '-dpng', ['-r' num2str(500)]); close(f1);

    % second plot: kkt compare legend
    rows = 50; cols = 500; linewidth = 3;
    LineColors = linspecer(2, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f3 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    legend_1 = "$r_{\max}$: Baseline"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$r_{\max}$: Penalty Wave"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend({legend_1, legend_2}, 'Interpreter', 'latex', 'NumColumns', 2, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f3, img_prefix + "sigwave_kkt_comp_legend.png", '-dpng', ['-r' num2str(500)]); close(f3);
end
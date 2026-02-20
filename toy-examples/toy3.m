clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

addpath(genpath("../utils/"));
addpath(genpath("../admm-single-light/"));

%% switches
if_taskI = false;
if_taskIV = false;
if_taskVI = true;
if_taskVII = false;

%% SDP data
SDP_data = generate_SDP_data();

%% task I: sweep t
if if_taskI 
    sig = 1;
    h = 1;
    eps = 0; % eps >= 0
    t_list = 10.^([-1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5]);
    % t_list = 10.^([-4]);
    maxiter = 1000;
    tol = 1e-14;
    data = cell(1, length(t_list));
    for ii = 1: length(t_list)
        t = t_list(ii);
        [Z_init, mdv_mats] = generate_init_guess(sig, h, eps, t);
        sig_data = set_sig(sig);
        iter_data = set_iter(Z_init, sig_data);
        iter_data.mdv_mats = mdv_mats;
        iter_data.t = t;
        [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data);
        data{ii}.t = t;
        data{ii}.Z_diff_norm_list = draw_data.Z_diff_norm_list;
        data{ii}.Z_diff_ang_list = draw_data.Z_diff_ang_list;
        data{ii}.mdvZ_relerr_list = draw_data.mdvZ_relerr_list;
    end
    draw_taskI(data, "./figs/toy3/");
end

%% task IV: update sigma
if if_taskIV
    maxiter = 1000;
    tol = -1;
    t = 1e-5;
    Ampl = 10;
    sig_list = Ampl.^linspace(0, 1, maxiter);
    h = 1; 
    eps = 0;
    [Z_init, mdv_mats] = generate_init_guess(sig_list(1), h, eps, t);
    sig_data = set_sig(sig_list(1));
    iter_data = set_iter(Z_init, sig_data);
    iter_data.mdv_mats = mdv_mats;
    iter_data.t = t;
    [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data, sig_list);
    
    % first plot: pinf, dinf, deltaX, deltaS
    rows = 300; cols = 300; linewidth = 1.5; ms = 12; xls = 15; yls = 15; fontsize = 15;
    LineColors = linspecer(4, 'qualitative'); LineColors = flipud(LineColors); 
    % Prepare figure
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    hold on; grid on;
    plot(log10(draw_data.sig_list), log10(draw_data.X_diff_norm_list), 'LineWidth', linewidth, 'Color', LineColors(1, :));
    loglog(log10(draw_data.sig_list), log10(draw_data.S_diff_norm_list), 'LineWidth', linewidth, 'Color', LineColors(2, :));
    loglog(log10(draw_data.sig_list), log10(draw_data.pinf_list), 'LineWidth', linewidth, 'Color', LineColors(3, :));
    loglog(log10(draw_data.sig_list), log10(draw_data.dinf_list), 'LineWidth', linewidth, 'Color', LineColors(4, :));
    xlabel('$\log_{10}(\sigma)$', 'FontSize', fontsize, 'Interpreter', 'latex');
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f1, sprintf("./figs/toy3/taskIV_tinv=%d.png", ceil(1/t)), '-dpng', ['-r' num2str(500)]); close(f1)
    % second plot: legend
    % second plot: original legend
    rows = 50; cols = 500; linewidth = 3;
    LineColors = linspecer(4, 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f2 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    legend_1 = "$\log_{10}(\parallel \Delta X^{(k)} \parallel_F)$"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$\log_{10}(\parallel \Delta S^{(k)} \parallel_F)$"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend_3 = "$\log_{10}(r_p^{(k)})$"; plot(nan, 'Color', LineColors(3, :), 'LineWidth', linewidth); 
    legend_4 = "$\log_{10}(r_d^{(k)})$"; plot(nan, 'Color', LineColors(4, :), 'LineWidth', linewidth); 
    legend({legend_1, legend_2, legend_3, legend_4}, 'Interpreter', 'latex', 'NumColumns', 4, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f2, "./figs/toy3/taskIV_legend.png", '-dpng', ['-r' num2str(500)]); close(f2);
end

%% task VI: update sig to escape the slow convergence region
if if_taskVI
    sig_change_list = [1e-2, 1e-1, 1, 1e1, 1e2];
    t = 1e-4;
    h = 1; 
    eps = 0;
    maxiter = 2000;
    tol = 1e-14;

    rows = 350; cols = 350; linewidth = 1.5; ms = 12; xls = 15; yls = 15; zls = 15; fontsize = 15;
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    grid on; hold on;
    LineColors = linspecer(length(sig_change_list), 'qualitative'); LineColors = flipud(LineColors); 
    
    for i = 1: length(sig_change_list)
        sig_list = [ones(1, maxiter/2), sig_change_list(i) * ones(1, maxiter/2)];
        sig = sig_list(1); sig_data = set_sig(sig);
        [Z_init, mdv_mats] = generate_init_guess(sig, h, eps, t);
        iter_data = set_iter(Z_init, sig_data);
        iter_data.mdv_mats = mdv_mats;
        iter_data.t = t;
        [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data, sig_list);
        kkt_list = max([draw_data.pinf_list; draw_data.dinf_list; draw_data.relgap_list]);
        kkt_list = min(1, kkt_list);
        plot(log10(kkt_list), 'LineWidth', linewidth, 'Color', LineColors(i, :));
    end

    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    ylabel('$\log_{10}(r^{(k)}_{\max})$', 'FontSize', fontsize, 'Interpreter', 'latex');
    axis tight;
    print(f1, "./figs/toy3/taskVI.png", '-dpng', ['-r' num2str(500)]); close(f1);

    % labels
    rows = 50; cols = 500; linewidth = 3;
    % Prepare figure
    f2 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;
    for i = 1: length(sig_change_list)
        n = round(log10(sig_change_list(i)));                 
        s = sprintf('$\\sigma = 10^{%d}$', n); 
        plot(nan, 'Color', LineColors(i, :), 'LineWidth', linewidth, 'DisplayName', s); 
    end
    legend('show', 'Interpreter', 'latex', 'NumColumns', length(sig_change_list), 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f2, "./figs/toy3/taskVI_legend.png", '-dpng', ['-r' num2str(500)]); close(f2);
end

%% task VII: spike
if if_taskVII
    t = 1e-4;
    h_list = [1.60, 1.50, 1.40, 1.30];
    eps_list = [1e-2, 1e-3, 1e-4, 1e-5];
    maxiter = 1000;
    tol = -1;
    
    rmdir("./figs/toy3/taskVII", 's');
    mkdir("./figs/toy3/taskVII");
    for i = 1: length(eps_list)
        for j = 1: length(h_list)
            eps = eps_list(i);
            h = h_list(j);
            sig_list = [ones(1, maxiter/2), 1 * ones(1, maxiter/2)];
            sig = sig_list(1); sig_data = set_sig(sig);
            [Z_init, mdv_mats] = generate_init_guess(sig, h, eps, t);
            iter_data = set_iter(Z_init, sig_data);
            iter_data.mdv_mats = mdv_mats;
            iter_data.t = t;
            [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data, sig_list);
            
            % kkt
            rows = 150; cols = 300; linewidth = 1.5; ms = 12; xls = 15; yls = 15; zls = 15; fontsize = 15;
            f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
            grid on; hold on;
            LineColors = linspecer(3, 'qualitative'); LineColors = flipud(LineColors);
            
            Znorm_diffdiff = abs(diff(draw_data.Z_diff_norm_list));
            plot(log10(Znorm_diffdiff), 'LineWidth', linewidth, 'Color', LineColors(1, :));
            kkt_list = max([draw_data.pinf_list; draw_data.dinf_list; draw_data.relgap_list]);
            kkt_list_diff = abs(diff(kkt_list));
            plot(log10(kkt_list_diff), 'LineWidth', linewidth, 'Color', LineColors(2, :));
            
            % title
            xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
            % Axes formatting
            ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
            ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
            axis tight;
            % Save & close
            print(f1, sprintf("./figs/toy3/taskVII/%d_%d_kkt.png", i, j), '-dpng', ['-r' num2str(500)]); close(f1);
            
            % angle
            title_size = 18;
            rows = 150; cols = 300; linewidth = 1.5; ms = 12; xls = 15; yls = 15; zls = 15; fontsize = 15;
            f2 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
            grid on; hold on;
            plot(log10(draw_data.Z_diff_ang_list), 'LineWidth', linewidth, 'Color', LineColors(3, :));
            n = round(log10(eps));
            title_name = sprintf("$\\epsilon = 10^{%d}, h = %.2f$", n, h);
            title(title_name, 'FontSize', title_size, 'Interpreter', 'latex');
            % Axes formatting
            ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
            ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
            axis tight;
            % Save & close
            print(f2, sprintf("./figs/toy3/taskVII/%d_%d_ang.png", i, j), '-dpng', ['-r' num2str(500)]); close(f2);
        end
    end

    % labels
    rows = 50; cols = 700; linewidth = 3;
    % Prepare figure
    f3 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); hold on;               
    s = "$\log_{10}(| \parallel \Delta Z^{(k+1)} \parallel_F - \parallel \Delta Z^{(k)} \parallel_F |)$";
    plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth, 'DisplayName', s); 
    s = "$\log_{10}(|r_{\max}^{(k+1)} - r_{\max}^{(k)}|)$";
    plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth, 'DisplayName', s); 
    s = "$\log_{10}(\angle (\Delta Z^{(k)}, \Delta Z^{(k+1)}) )$";
    plot(nan, 'Color', LineColors(3, :), 'LineWidth', linewidth, 'DisplayName', s); 
    legend('show', 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off'); set(gca, 'LooseInset', [0,0,0,0]); set(gca, 'Box', 'off'); set(gcf, 'Color', 'w'); axis off; 
    % Save & close
    print(f3, "./figs/toy3/taskVII_legend.png", '-dpng', ['-r' num2str(500)]); close(f3);
end



%% helper functions
function [Z_init, mats] = generate_init_guess(sig, h, eps, t)
    % Z_bar
    X_bar_mat = zeros(6); X_bar_mat(1,1) = 6;
    S_bar_mat = zeros(6); S_bar_mat(6,6) = 3;
    Z_bar_mat = X_bar_mat - sig * S_bar_mat;
    
    % H_bar
    a = 1;
    H11 = -a;
    H13 = -a/(2*sqrt(2));
    H22 = a;
    H14 = 1;
    H15 = h;
    HX_bar_mat = [
        H11, 0  , H13, H14, H15, 0;
        0  , H22, 0  , 0  , 0  , 0;
        H13, 0  , 0  , 0  , 0  , 0;
        H14, 0  , 0  , 0  , 0  , 0;
        H15, 0  , 0  , 0  , 0  , 0;
        0  , 0  , 0  , 0  , 0  , 0;
    ];
    H26 = 1;
    H36 = 1;
    H46 = 1;
    H56 = 1;
    HS_bar_mat = [
        0, 0  , 0  , 0  , 0  , 0  ;
        0, 0  , 0  , 0  , 0  , H26;
        0, 0  , 0  , 0  , 0  , H36;
        0, 0  , 0  , -eps  , 0  , H46;
        0, 0  , 0  , 0  , -1 , H56;
        0, H26, H36, H46, H56, 1+eps  ;
    ];
    H_bar_mat = HX_bar_mat + HS_bar_mat;
    
    Z_init = svec_single(Z_bar_mat + t * H_bar_mat);
    mats.X_bar = svec_single(X_bar_mat); mats.X_bar_mat = X_bar_mat;
    mats.S_bar = svec_single(S_bar_mat); mats.S_bar_mat = S_bar_mat;
    mats.Z_bar = svec_single(Z_bar_mat); mats.Z_bar_mat = Z_bar_mat;
    mats.HX_bar = svec_single(HX_bar_mat); mats.HX_bar_mat = HX_bar_mat;
    mats.HS_bar = svec_single(HS_bar_mat); mats.HS_bar_mat = HS_bar_mat;
    mats.H_bar = svec_single(H_bar_mat); mats.H_bar_mat = H_bar_mat;

    mats.X = [
        -4/(9*sig)        , -2*sqrt(2)/(9*sig),         0, -eps/3, -h/3, 0;
        -2*sqrt(2)/(9*sig),          2/(9*sig), 4/(9*sig),      0,    0, 0;
                         0,          4/(9*sig), 2/(9*sig),      0,    0, 0;
                    -eps/3,                  0,         0,      0,    0, 0;
                      -h/3,                  0,         0,      0,    0, 0;
                         0,                  0,         0,      0,    0, 0;
    ];
    if abs(h) <= sqrt(2) || (abs(h) > sqrt(2) && eps > 0)
        mats.S = -1/sig * [ 
            0,  0,  0,             0,                 0,                0;
            0,  0,  0,             0,                 0,     -2/(3*sig);
            0,  0,  0,             0,     sqrt(2)*h/12,                0;
            0,  0,  0,   (h^2 - 2)/9,            -h/3,   2*eps/(3*sig);
            0,  0,  sqrt(2)*h/12, -h/3,  -(2*h^2 - 1)/9,               0;
            0, -2/(3*sig), 0,  2*eps/(3*sig),     0,      (h^2 + 1)/9;
        ];
    else
        mats.S = -1/sig * [ 
            0,  0,  0,             0,                 0,                0;
            0,  0,  0,             0,                 0,     -2/(3*sig);
            0,  0,  0,             0,     sqrt(2)*h/12,                0;
            0,  0,  0,   0,            -h/3,   2*eps/(3*sig);
            0,  0,  sqrt(2)*h/12, -h/3,  -h^2/6,               0;
            0, -2/(3*sig), 0,  2*eps/(3*sig),     0,      h^2/6;
        ];
    end
    mats.Z = mats.X - sig * mats.S;
end 

function SDP_data = generate_SDP_data()
    % C
    C = zeros(6); C(4,4) = 1; C(5,5) = 1; C(6,6) = 1; C = svec_single(C);
    
    % A
    q1 = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
    q2 = [1/sqrt(2); -1/sqrt(2); 0];
    q3 = [1/sqrt(6); 1/sqrt(6); -2/sqrt(6)];
    Q = [q1, q2, q3];
    Amats_cellarr = cell(0);
    A = eye(6); Amats_cellarr{end+1} = A;
    A = [Q'*diag([1,-1,0])*Q, zeros(3); zeros(3), zeros(3)]; Amats_cellarr{end+1} = A;
    A = [Q'*diag([1,0,-1])*Q, zeros(3); zeros(3), zeros(3)]; Amats_cellarr{end+1} = A;
    A = zeros(6); A(4,4) = 1; A(5,5) = -1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(5,5) = 1; A(6,6) = -1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(2,4) = 1; A(4,2) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(2,5) = 1; A(5,2) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(2,6) = 1; A(6,2) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(3,4) = 1; A(4,3) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(3,5) = 1; A(5,3) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(3,6) = 1; A(6,3) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(4,5) = 1; A(5,4) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(4,6) = 1; A(6,4) = 1; Amats_cellarr{end+1} = A;
    A = zeros(6); A(5,6) = 1; A(6,5) = 1; Amats_cellarr{end+1} = A;
    % test: (+, -) block's constraint effect
    A = zeros(6); A(1,6) = 1; A(6,1) = 1; Amats_cellarr{end+1} = A;
    At = assmeble_Amats(Amats_cellarr);

    % b
    b = zeros(length(Amats_cellarr), 1);
    b(1) = 6;

    SDP_data.C = C;
    SDP_data.At = At;
    SDP_data.A = At';
    SDP_data.b = b;
    SDP_data.PA = SDP_data.At * inv(SDP_data.A * SDP_data.At) * SDP_data.A;
    SDP_data.PAp = eye(21) - SDP_data.PA;
    SDP_data.xi = SDP_data.A \ SDP_data.b;
end

function At = assmeble_Amats(Amats_cellarr) 
    At = svec_single(Amats_cellarr{1});
    for i = 2: length(Amats_cellarr)
        At = [At, svec_single(Amats_cellarr{i})];
    end
end

function [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data, sig_list)
    t = iter_data.t;
    mdv_mats_theory = iter_data.mdv_mats;

    if nargin == 6
        sig_data.sig = sig_list(1);
    end

    [iter_data_prev, SDP_data] = scale_data(iter_data, sig_data, SDP_data, false); 
    iter_data = three_step_ADMM(iter_data_prev, sig_data, SDP_data); 
    iter_data_next = three_step_ADMM(iter_data, sig_data, SDP_data); 
    log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
    draw_data = set_draw(0, [], sig_data, log_data, true);

    mats = struct();
    mats.X = cell(1, maxiter);
    mats.S = cell(1, maxiter);
    mats.Z = cell(1, maxiter);

    diff_mats = struct();
    diff_mats.X = cell(1, maxiter);
    diff_mats.S = cell(1, maxiter);
    diff_mats.Z = cell(1, maxiter);

    draw_data.mdvZ_relerr_list = [];
    draw_data.mdvX_relerr_list = [];
    draw_data.mdvS_relerr_list = [];

    for iter = 1: maxiter
        iter_data_prev = iter_data;
        iter_data = iter_data_next;
        iter_data_next = three_step_ADMM(iter_data_next, sig_data, SDP_data);
        log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
        draw_data = set_draw(iter, draw_data, sig_data, log_data);
        
        mats.X{iter} = smat_single(iter_data_next.X);
        mats.S{iter} = smat_single(iter_data_next.S);
        mats.Z{iter} = smat_single(iter_data_next.Z);
        Z_diff_mat = smat_single(iter_data_next.Z - iter_data.Z);
        X_diff_mat = smat_single(iter_data_next.X - iter_data.X);
        S_diff_mat = smat_single(iter_data_next.S - iter_data.S);
        diff_mats.X{iter} = X_diff_mat;
        diff_mats.S{iter} = S_diff_mat;
        diff_mats.Z{iter} = Z_diff_mat;

        % disp(2 / (t^2) * Z_diff_mat);

        % calculate relative error between theory and practice
        mdvZ_relerr = norm(mdv_mats_theory.Z - Z_diff_mat/(0.5*t^2), 'fro') / norm(Z_diff_mat/(0.5*t^2), 'fro');
        mdvX_relerr = norm(mdv_mats_theory.X - X_diff_mat/(0.5*t^2), 'fro') / norm(X_diff_mat/(0.5*t^2), 'fro');
        mdvS_relerr = norm(mdv_mats_theory.S - S_diff_mat/(0.5*t^2), 'fro') / norm(S_diff_mat/(0.5*t^2), 'fro');
        draw_data.mdvZ_relerr_list = [draw_data.mdvZ_relerr_list, mdvZ_relerr];
        draw_data.mdvX_relerr_list = [draw_data.mdvX_relerr_list, mdvX_relerr];
        draw_data.mdvS_relerr_list = [draw_data.mdvS_relerr_list, mdvS_relerr];

        if log_data.kkt < tol 
            break;
        end

        if nargin == 6
            sig_data.sig = sig_list(iter);
        end

        if mod(iter, 1) == 0 
            fprintf("iter: %d, [%3.2e, %3.2e, %d], [%3.2e | %3.2e | %3.2e], [%3.2e, %3.2e], %3.2e \n",...
                        iter, log_data.Z_diff_norm, log_data.Z_diff_ang, log_data.Z_rank, ...
                        log_data.pinf, log_data.dinf, log_data.relgap, ...
                        log_data.pobj, log_data.dobj, sig_data.sig);
        end
    end
end

% symbolic helper
function symbolic_helper()
    % a represents sigma
    syms a b real;
    assumeAlso(a > 0);
    sqrt2 = sqrt(sym(2));
    sqrt3 = sqrt(sym(3));
    sqrt6 = sqrt(sym(6));
    sig = a;
    sigp = b;

    % define Z's eigenvalues
    % there is only one positive eigenvalue and one negative eigenvalue
    mu_p = 6;
    mu_n = -3*sig * sigp/sig;

    % define H
    coeff = 1;
    H11 = -coeff;
    H13 = -coeff/(2*sqrt2);
    H22 = coeff;
    H14 = 1;
    H15 = 1;
    HX_bar_mat = [
        H11, 0  , H13, H14, H15, 0;
        0  , H22, 0  , 0  , 0  , 0;
        H13, 0  , 0  , 0  , 0  , 0;
        H14, 0  , 0  , 0  , 0  , 0;
        H15, 0  , 0  , 0  , 0  , 0;
        0  , 0  , 0  , 0  , 0  , 0;
    ];
    H26 = 1;
    H36 = 1;
    H46 = 1;
    H56 = 1;
    HS_bar_mat = [
        0, 0  , 0  , 0  , 0  , 0  ;
        0, 0  , 0  , 0  , 0  , H26;
        0, 0  , 0  , 0  , 0  , H36;
        0, 0  , 0  , 0  , 0  , H46;
        0, 0  , 0  , 0  , -1 , H56;
        0, H26, H36, H46, H56, 1  ;
    ];
    H = HX_bar_mat + sigp / sig * HS_bar_mat;

    % define E and Ep
    H_pp = H(1:1, 1:1); H_p0 = H(1:1, 2:5); H_pn = H(1:1, 6:6);
    H_0p = H(2:5, 1:1); H_00 = H(2:5, 2:5); H_0n = H(2:5, 6:6);
    H_np = H(6:6, 1:1); H_n0 = H(6:6, 2:5); H_nn = H(6:6, 6:6);
    H_00_p = H_00; H_00_p(4,4) = 0;
    H_00_n = H_00; H_00_n(1,1) = 0;
    E_pp = 0;
    E_p0 = -2 * 1/mu_p * H_p0 * -H_00_n;
    E_pn = 2/(mu_p - mu_n) * H_p0 * H_0n;
    E_0p = E_p0';
    E_00 = 2 * 1/mu_p * H_0p * H_p0;
    E_0n = 2/(-mu_n) * H_00_p * H_0n;
    E_np = E_pn';
    E_n0 = E_0n';
    E_nn = 0;
    E = [
        E_pp, E_p0, E_pn;
        E_0p, E_00, E_0n;
        E_np, E_n0, E_nn;
    ];
    Ep_00 = 2 * 1/mu_n * H_0n * H_n0;
    Ep = [
        -E_pp, -E_p0, -E_pn;
        -E_0p, Ep_00, -E_0n;
        -E_np, -E_n0, -E_nn;
    ];

    E = simplify(E);
    Ep = simplify(Ep);

    % further computation
    q1 = [1/sqrt3; 1/sqrt3; 1/sqrt3];
    q2 = [1/sqrt2; -1/sqrt2; 0];
    q3 = [1/sqrt6; 1/sqrt6; -2/sqrt6];
    Q = [q1, q2, q3];
    Tmp = Q * [
        0, 0, 0;
        0, -2/(3*sig), -2/(3*sig);
        0, -2/(3*sig), -2/(3*sig); 
    ] * Q';
    Tmp = simplify(expand(Tmp), 'Steps', 100);
    Tmp(1,1) = 0;
    Tmp(2,2) = 0;
    Tmp(3,3) = 0;
    Tmp = Q' * Tmp * Q;
    Tmp = simplify(Tmp, 'Steps', 100);

    Tmp1 = Q * [
        0, 0, 0;
        0, 1, 1;
        0, 1, 1; 
    ] * Q';
    Tmp1 = simplify(Tmp1, 'Steps', 100);
    Tmp1(1,1) = 0;
    Tmp1(2,2) = 0;
    Tmp1(3,3) = 0;
    Tmp1 = Q' * Tmp1 * Q;
    Tmp1 = simplify(Tmp1, 'Steps', 100);
end

% symbolic helper
function Z_init = symbolic_helper_1(t_num, sig_num, h_num, eps_num)
    % a represents sigma
    syms a b real;
    % h is H15, e is epsilon
    syms h e real;
    assumeAlso(a > 0);
    assumeAlso(e >= 0);
    sqrt2 = sqrt(sym(2));
    sqrt3 = sqrt(sym(3));
    sqrt6 = sqrt(sym(6));
    sig = a;
    sigp = a;
    eps = e;

    % A
    q1 = [1/sqrt3; 1/sqrt3; 1/sqrt3];
    q2 = [1/sqrt2; -1/sqrt2; 0];
    q3 = [1/sqrt6; 1/sqrt6; -2/sqrt6];
    Q = [q1, q2, q3];
    Amats_cellarr = cell(0);
    A = eye(6); Amats_cellarr{end+1} = A;
    A = [Q'*diag([1,-1,0])*Q, zeros(3); zeros(3), zeros(3)]; Amats_cellarr{end+1} = sym(A);
    A = [Q'*diag([1,0,-1])*Q, zeros(3); zeros(3), zeros(3)]; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(4,4) = 1; A(5,5) = -1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(5,5) = 1; A(6,6) = -1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(2,4) = 1; A(4,2) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(2,5) = 1; A(5,2) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(2,6) = 1; A(6,2) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(3,4) = 1; A(4,3) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(3,5) = 1; A(5,3) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(3,6) = 1; A(6,3) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(4,5) = 1; A(5,4) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(4,6) = 1; A(6,4) = 1; Amats_cellarr{end+1} = sym(A);
    A = zeros(6); A(5,6) = 1; A(6,5) = 1; Amats_cellarr{end+1} = sym(A);

    % define Z's eigenvalues
    % there is only one positive eigenvalue and one negative eigenvalue
    mu_p = 6;
    mu_n = -3*sig * sigp/sig;

    % define H
    coeff = 1;
    H11 = -coeff;
    H13 = -coeff/(2*sqrt2);
    H22 = coeff;
    H14 = 1;
    H15 = h;
    HX_bar_mat = [
        H11, 0  , H13, H14, H15, 0;
        0  , H22, 0  , 0  , 0  , 0;
        H13, 0  , 0  , 0  , 0  , 0;
        H14, 0  , 0  , 0  , 0  , 0;
        H15, 0  , 0  , 0  , 0  , 0;
        0  , 0  , 0  , 0  , 0  , 0;
    ];
    H26 = 1;
    H36 = 1;
    H46 = 1;
    H56 = 1;
    HS_bar_mat = [
        0, 0  , 0  , 0   , 0  ,   0  ;
        0, 0  , 0  , 0   , 0  ,   H26;
        0, 0  , 0  , 0   , 0  ,   H36;
        0, 0  , 0  , -eps, 0  ,   H46;
        0, 0  , 0  , 0   , -1 ,   H56;
        0, H26, H36, H46 , H56, 1+eps;
    ];
    H = HX_bar_mat + sigp / sig * HS_bar_mat;

    % define E and Ep
    H_pp = H(1:1, 1:1); H_p0 = H(1:1, 2:5); H_pn = H(1:1, 6:6);
    H_0p = H(2:5, 1:1); H_00 = H(2:5, 2:5); H_0n = H(2:5, 6:6);
    H_np = H(6:6, 1:1); H_n0 = H(6:6, 2:5); H_nn = H(6:6, 6:6);
    H_00_p = H_00; H_00_p(4,4) = 0;
    H_00_n = H_00; H_00_n(1,1) = 0;
    E_pp = 0;
    E_p0 = -2 * 1/mu_p * H_p0 * -H_00_n;
    E_pn = 2/(mu_p - mu_n) * H_p0 * H_0n;
    E_0p = E_p0';
    E_00 = 2 * 1/mu_p * H_0p * H_p0;
    E_0n = 2/(-mu_n) * H_00_p * H_0n;
    E_np = E_pn';
    E_n0 = E_0n';
    E_nn = 0;
    E = [
        E_pp, E_p0, E_pn;
        E_0p, E_00, E_0n;
        E_np, E_n0, E_nn;
    ];
    Ep_00 = 2 * 1/mu_n * H_0n * H_n0;
    Ep = [
        -E_pp, -E_p0, -E_pn;
        -E_0p, Ep_00, -E_0n;
        -E_np, -E_n0, -E_nn;
    ];

    E = simplify(E);
    Ep = simplify(Ep);
    
    disp(latex(H));
    disp(latex(E));
    disp(latex(Ep));

    % return Z_init
    H_num = double(subs(H, [a, b, h, e], [sig_num, sig_num, h_num, eps_num]));
    Zbar_num = sym(zeros(6));
    Zbar_num(1,1) = mu_p; Zbar_num(6, 6) = mu_n;
    Zbar_num = double(subs(Zbar_num, [a, b, h, e], [sig_num, sig_num, h_num, eps_num]));
    Z_init = svec_single(Zbar_num + t_num * H_num);
end

function b = AX(Amats_cellarr, Xmat)
    b = zeros(length(Amats_cellarr), 1);
    for i = 1: length(b)
        b(i) = trace(Amats_cellarr{i} * Xmat);
    end
end

function draw_taskI(data, img_prefix)
    % first plot: Znorm compare
    rows = 450; cols = 450; linewidth = 1.5; ms = 12; xls = 15; yls = 15; fontsize = 15;
    LineColors = linspecer(length(data), 'qualitative'); LineColors = flipud(LineColors);
    % Prepare figure
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % plot
    hold on; grid on; 
    for ii = 1: length(data)
        legend_name = sprintf("$\\log_{10}(t) = %.1f$", log10(data{ii}.t)); 
        plot(log10(data{ii}.Z_diff_norm_list), 'LineWidth', linewidth, 'Color', LineColors(ii, :), 'DisplayName', legend_name);
    end
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    ylabel('$\log_{10}(\parallel \Delta Z^{(k)} \parallel_F)$', 'FontSize', fontsize, 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex', 'NumColumns', ceil(length(data)/3), 'Location', 'northoutside', 'FontSize', fontsize);
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f1, img_prefix + "taskI_Znorm_comp.png", '-dpng', ['-r' num2str(500)]); close(f1);

    % second plot: mdvZ_relerr compare
    rows = 450; cols = 450; linewidth = 1.5;
    % Prepare figure
    f2 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % plot
    hold on; grid on; 
    for ii = 1: length(data)
        legend_name = sprintf("$\\log_{10}(t) = %.1f$", log10(data{ii}.t)); 
        plot(log10(min(data{ii}.mdvZ_relerr_list, 1)), 'LineWidth', linewidth, 'Color', LineColors(ii, :), 'DisplayName', legend_name);
    end
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    ylabel('$\log_{10}(\parallel \frac{1}{2} \phi(\bar{Z}; Z^{(0)} - \bar{Z}) - \Delta Z^{(k)} \parallel_F / \parallel \Delta Z^{(k)} \parallel_F)$', 'FontSize', fontsize, 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex', 'NumColumns', ceil(length(data)/3), 'Location', 'northoutside', 'FontSize', fontsize);
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f2, img_prefix + "taskI_mdvZ_relerr_comp.png", '-dpng', ['-r' num2str(500)]); close(f2);

    % third plot: Z_diff_ang compare
    rows = 450; cols = 450; linewidth = 1.5;
    % Prepare figure
    f4 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % plot
    hold on; grid on; 
    for ii = 1: length(data)
        legend_name = sprintf("$\\log_{10}(t) = %.1f$", log10(data{ii}.t)); 
        plot(log10(data{ii}.Z_diff_ang_list), 'LineWidth', linewidth, 'Color', LineColors(ii, :), 'DisplayName', legend_name);
    end
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    ylabel('$\log_{10}(\angle (\Delta Z^{(k)}, \Delta Z^{(k+1)}) )$', 'FontSize', fontsize, 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex', 'NumColumns', ceil(length(data)/3), 'Location', 'northoutside', 'FontSize', fontsize);
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f4, img_prefix + "taskI_Zang_comp.png", '-dpng', ['-r' num2str(500)]); close(f4);
end








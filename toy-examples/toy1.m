clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

addpath(genpath("../utils/"));
addpath(genpath("../admm-single-light/"));

%% switches
if_taskI = false;
if_taskII = false;
if_taskIII = false;
if_taskIV = false;
if_taskV = true;
if_taskVI = false;

%% SDP data
SDP_data = generate_SDP_data();

%% task I: sweep t
if if_taskI 
    sig = 1;
    a = 1; % a >= 0
    b = 1;
    t_list = 10.^([-1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5]);
    % t_list = 10.^([-4]);
    maxiter = 1000;
    tol = 1e-14;
    data = cell(1, length(t_list));
    for ii = 1: length(t_list)
        t = t_list(ii);
        [Z_init, mdv_mats] = generate_init_guess(sig, a, b, t);
        sig_data = set_sig(sig);
        iter_data = set_iter(Z_init, sig_data);
        iter_data.mdv_mats = mdv_mats;
        iter_data.t = t;
        [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data);
        data{ii}.t = t;
        data{ii}.Z_diff_norm_list = draw_data.Z_diff_norm_list;
        data{ii}.Z_diff_ang_list = draw_data.Z_diff_ang_list;
        data{ii}.mdvZ_relerr_list = draw_data.mdvZ_relerr_list;
        data{ii}.mdvZ_instant_relerr_list = draw_data.mdvZ_instant_relerr_list;
    end
    draw_taskI(data, "./figs/toy1/");
end

%% task II: field plots in 3D
if if_taskII
    t_list = 10.^([-1, -2, -3, -4]);
    sig = 1;
    sample_num = 10;
    % sample a^2 + b^2 = 1 with a >= 0
    th_list = linspace(-pi/2, pi/2, sample_num);
    t = t_list(4);
    figure;
    for i = 1: length(th_list)
        th = th_list(i);
        a = cos(th); b = sin(th);
        [Z_init, mdv_mats] = generate_init_guess(sig, a, b, t);
        sig_data = set_sig(sig);
        iter_data = set_iter(Z_init, sig_data);
        iter_data.mdv_mats = mdv_mats;
        iter_data.t = t;
        [draw_data, mats, diff_mats] = run_ADMM(50, -1, iter_data, sig_data, SDP_data);
        [x_list, y_list, z_list] = get_Z_3D_lists(mats);
        % set color
        % tmp = linspace(0,1,length(x_list))';
        % rgb = [tmp, zeros(length(x_list),1), 1-tmp];
        handle = scatter3(x_list(:), y_list(:), z_list(:), 20, 'filled', 'MarkerEdgeColor','none'); hold on; axis equal;
        handle.AlphaData = linspace(1,0,length(x_list)).'; handle.AlphaDataMapping = 'none'; handle.MarkerFaceAlpha = 'flat'; 
    end
end

%% task III: theoretical plots
if if_taskIII
    sig = 1;
    t = 1e-2;
    N_a = 10;
    N_b = 20;
    a_list = linspace(0, 1, N_a);
    b_list = linspace(-1, 1, N_b);
    [x, y, z] = ndgrid(t*a_list, t*b_list, t*b_list-sig);
    u = zeros(size(x)); v = zeros(size(x)); w = zeros(size(x));
    for i = 1: N_a 
        for j = 1: N_b
            k = N_b - j + 1;
            a = a_list(i); b = b_list(j);
            [Z_init, mdv_mats] = generate_init_guess(sig, a, b, t);
            Z_mat = [0, 0; 0, -sig] + t * [a, b; b , -b];
            Z_diff_mat = 0.5 * t^2 * mdv_mats.Z;
            u(i, j, k) = Z_diff_mat(1,1);
            v(i, j, k) = Z_diff_mat(1,2);
            w(i, j, k) = Z_diff_mat(2,2);
        end
    end
    
    % plot
    rows = 500; cols = 500; linewidth = 1.5; ms = 12; xls = 15; yls = 15; zls = 15; fontsize = 15;
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on');
    % f1 = figure;
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % vector field
    hq = quiver3(x,y,z,u,v,w, 'Color', 'r', 'AutoScale','on', 'AutoScaleFactor', 1.1); axis equal; grid on; hold on;
    hq.LineWidth = 1.0; hq.MaxHeadSize = 0.6; axis equal;
    hq.DisplayName = "$\phi_{Z}(\bar{Z}, \cdot)$";
    % ker delta'_Z
    na = 10; nb = 10;
    [A, B] = meshgrid(linspace(0,t,na), linspace(-t,t,nb));  
    Xp = A; Yp = B; Zp = -sig-B; hp = surf(Xp, Yp, Zp);
    set(hp, 'FaceColor','b', 'FaceAlpha', 0.12, 'EdgeColor', 'none'); % small alpha
    hp.DisplayName = "$\mathcal{C}(\bar{Z})$";
    % optimal set
    hOpt = plot3([0, 1.2*t], [0, 0], [-sig, -sig], 'LineWidth', 1.5, 'Color', [0, 0, 0]);
    hOpt.DisplayName = "$\mathcal{Z}_\star$";
    % Z_bar 
    hZbar = scatter3(0, 0, -sig, 50, 'black', 'filled');
    hZbar.DisplayName = "$\bar{Z}$";
    % legend
    legend([hq, hp, hOpt, hZbar], 'Location', 'best', 'Interpreter', 'latex', 'FontSize', fontsize); 
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); 
    ax.XAxis.FontSize = xls; ax.YAxis.FontSize = yls; ax.ZAxis.FontSize = zls;
    ax.GridLineWidth  = 0.4; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    % view
    view(-45, 60);
    print(f1, "./figs/toy1/taskIII_theory.png", '-dpng', ['-r' num2str(500)]); close(f1);
end

%% task IV: update sigma
if if_taskIV
    maxiter = 1000;
    tol = -1;
    t = 1e-5;
    Ampl = 10;
    sig_list = Ampl.^linspace(0, 1, maxiter);
    a = 1; b = 1;
    [Z_init, mdv_mats] = generate_init_guess(sig_list(1), a, b, t);
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
    print(f1, sprintf("./figs/toy1/taskIV_tinv=%d.png", ceil(1/t)), '-dpng', ['-r' num2str(500)]); close(f1)
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
    print(f2, "./figs/toy1/taskIV_legend.png", '-dpng', ['-r' num2str(500)]); close(f2);
end


%% task V: start from outside
if if_taskV
    t = 1e-5;
    a = 0; 
    b = 0;
    maxiter = 1000;
    tol = -1;
    sig_list = [ones(1, maxiter/2), 1e0 * ones(1, maxiter/2)];
    sig = sig_list(1); sig_data = set_sig(sig);
    [Z_init, mdv_mats] = generate_init_guess(sig, a, b, t);

    % plot
    rows = 500; cols = 500; linewidth = 1.5; ms = 12; xls = 15; yls = 15; zls = 15; fontsize = 15;
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    axis equal; grid on; hold on; view(3);

    x_list = [-2, -1, 1, 2];
    y_list = [-2, -1, 1, 2];
    z_list = [-2, -1, 1, 2];
    
    for x = x_list
        for y = y_list
            for z = z_list
                iter_data = set_iter(Z_init, sig_data);
                H = [x, y; y, z];
                iter_data.noise = t * svec_single(H);
                iter_data.mdv_mats = mdv_mats;
                iter_data.t = t;
                [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data, sig_list);
            
                % vector field
                [Z11_list, Z12_list, Z22_list] = get_Z_3D_lists(mats);
                hq = plot3(Z11_list, Z12_list, Z22_list); 
                hq.LineWidth = 1.0; 
                hq.Color = 'r';
                hq.DisplayName = "$Z^{(k)}$";
                hs = scatter3(Z11_list, Z12_list, Z22_list, 10, 'r', 'filled', 'MarkerEdgeColor','none');
            end
        end
    end

    % ker delta'_Z
    na = 10; nb = 10;
    [A, B] = meshgrid(linspace(0,3*t,na), linspace(-3*t,3*t,nb));  
    Xp = A; Yp = B; Zp = -sig-B; hp = surf(Xp, Yp, Zp);
    set(hp, 'FaceColor','b', 'FaceAlpha', 0.12, 'EdgeColor', 'none'); % small alpha
    hp.DisplayName = "$\mathcal{C}(\bar{Z})$";
    % optimal set
    hOpt = plot3([0, 3.2*t], [0, 0], [-sig, -sig], 'LineWidth', 1.5, 'Color', [0, 0, 0]);
    hOpt.DisplayName = "$\mathcal{Z}_\star$";
    % Z_bar 
    hZbar = scatter3(0, 0, -sig, 50, 'black', 'filled');
    hZbar.DisplayName = "$\bar{Z}$";
    % legend
    legend([hq, hp, hOpt, hZbar], 'Location', 'best', 'Interpreter', 'latex', 'FontSize', fontsize); 
    % Axes formatting
    ax = gca; 
    ax.XAxis.FontSize = xls; ax.YAxis.FontSize = yls; ax.ZAxis.FontSize = zls;
    ax.GridLineWidth  = 0.4; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    ax.XLabel.String = '$Z_{11}$'; ax.XLabel.Interpreter = 'latex';
    ax.YLabel.String = '$Z_{12}$'; ax.YLabel.Interpreter = 'latex';
    ax.ZLabel.String = '$Z_{22}$'; ax.ZLabel.Interpreter = 'latex';
    set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); 
    % view
    view(-45, 60);
    % save
    filename = sprintf("./figs/toy1/taskV_tinv=%d.png", ceil(1/t));
    print(f1, filename, '-dpng', ['-r' num2str(500)]); close(f1);
end

%% task VI: update sig to escape the slow convergence region
if if_taskVI
    sig_change_list = [1e-2, 1e-1, 1, 1e1, 1e2];
    t = 1e-4;
    a = 1; 
    b = 1;
    maxiter = 2000;
    tol = 1e-14;

    rows = 350; cols = 350; linewidth = 1.5; ms = 12; xls = 15; yls = 15; zls = 15; fontsize = 15;
    f1 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    grid on; hold on;
    LineColors = linspecer(length(sig_change_list), 'qualitative'); LineColors = flipud(LineColors); 
    
    for i = 1: length(sig_change_list)
        sig_list = [ones(1, maxiter/2), sig_change_list(i) * ones(1, maxiter/2)];
        sig = sig_list(1); sig_data = set_sig(sig);
        [Z_init, mdv_mats] = generate_init_guess(sig, a, b, t);
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
    print(f1, "./figs/toy1/taskVI.png", '-dpng', ['-r' num2str(500)]); close(f1);

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
    print(f2, "./figs/toy1/taskVI_legend.png", '-dpng', ['-r' num2str(500)]); close(f2);
end





%% helper functions
function SDP_data = generate_SDP_data()
    % C
    C = [
        0, 0;
        0, 1;
    ];
    C = svec_single(C);
    
    % A
    A1 = [
        0, 1;
        1, -1;
    ];
    A1 = svec_single(A1);
    At = A1;

    % b
    b = 0;

    SDP_data.C = C;
    SDP_data.At = At;
    SDP_data.A = At';
    SDP_data.b = b;
end

function [Z_init, mdv_mats] = generate_init_guess(sig, a, b, t)
    Z_bar = [
        0, 0;
        0, -sig;
    ];
    H_bar = [
       a, b;
       b, -b;
    ];
    Z_init = svec_single(Z_bar + t * H_bar);

    mdv_mats.Z = [
        2*b^2/sig, -4/3*a*b/sig;
        -4/3*a*b/sig, 4/3*a*b/sig;
    ];
    mdv_mats.X = [
        2*b^2/sig, 0;
        0, 0;
    ];
    mdv_mats.S = [
        0,  4/3*a*b/(sig^2);
        4/3*a*b/(sig^2), -4/3*a*b/(sig^2);
    ];
end

function mdv_mats = instant_theory_mats(sig, Z)
    Z_bar = [
        0, 0;
        0, -sig;
    ];
    H = Z - Z_bar;
    x = H(1,1); y = H(1,2); z = H(2,2);
    a = x; b = 1/3*(2*y-z);

    mdv_mats.Z = [
        2*b^2/sig, -4/3*a*b/sig;
        -4/3*a*b/sig, 4/3*a*b/sig;
    ];
    mdv_mats.X = [
        2*b^2/sig, 0;
        0, 0;
    ];
    mdv_mats.S = [
        0,  4/3*a*b/(sig^2);
        4/3*a*b/(sig^2), -4/3*a*b/(sig^2);
    ];
end

function [draw_data, mats, diff_mats] = run_ADMM(maxiter, tol, iter_data, sig_data, SDP_data, sig_list)
    t = iter_data.t;
    mdv_mats_theory = iter_data.mdv_mats;
    iter_data_0 = iter_data;

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

    draw_data.mdvZ_instant_relerr_list = [];
    draw_data.mdvX_instant_relerr_list = [];
    draw_data.mdvS_instant_relerr_list = [];

    % inject noise
    if isfield(iter_data_0, 'noise')
        Z_new = iter_data_next.Z + iter_data_0.noise;
        iter_data_next = set_iter(Z_new, sig_data, SDP_data);
    end

    for iter = 1: maxiter
        mats.X{iter} = smat_single(iter_data_next.X);
        mats.S{iter} = smat_single(iter_data_next.S);
        mats.Z{iter} = smat_single(iter_data_next.Z);
        Z_diff_mat = smat_single(iter_data_next.Z - iter_data.Z);
        X_diff_mat = smat_single(iter_data_next.X - iter_data.X);
        S_diff_mat = smat_single(iter_data_next.S - iter_data.S);
        diff_mats.X{iter} = X_diff_mat;
        diff_mats.S{iter} = S_diff_mat;
        diff_mats.Z{iter} = Z_diff_mat;

        % calculate relative error between theory and practice
        mdvZ_relerr = norm(mdv_mats_theory.Z - Z_diff_mat/(0.5*t^2), 'fro') / norm(Z_diff_mat/(0.5*t^2), 'fro');
        mdvX_relerr = norm(mdv_mats_theory.X - X_diff_mat/(0.5*t^2), 'fro') / norm(X_diff_mat/(0.5*t^2), 'fro');
        mdvS_relerr = norm(mdv_mats_theory.S - S_diff_mat/(0.5*t^2), 'fro') / norm(S_diff_mat/(0.5*t^2), 'fro');
        draw_data.mdvZ_relerr_list = [draw_data.mdvZ_relerr_list, mdvZ_relerr];
        draw_data.mdvX_relerr_list = [draw_data.mdvX_relerr_list, mdvX_relerr];
        draw_data.mdvS_relerr_list = [draw_data.mdvS_relerr_list, mdvS_relerr];

        % calculate relative error between theory and practice
        mdv_mats_instant_theory = instant_theory_mats(sig_data.sig, mats.Z{iter});
        mdvZ_instant_relerr = norm(0.5*mdv_mats_instant_theory.Z - Z_diff_mat, 'fro') / norm(Z_diff_mat, 'fro');
        mdvX_instant_relerr = norm(0.5*mdv_mats_instant_theory.X - X_diff_mat, 'fro') / norm(X_diff_mat, 'fro');
        mdvS_instant_relerr = norm(0.5*mdv_mats_instant_theory.S - S_diff_mat, 'fro') / norm(S_diff_mat, 'fro');
        draw_data.mdvZ_instant_relerr_list = [draw_data.mdvZ_instant_relerr_list, mdvZ_instant_relerr];
        draw_data.mdvX_instant_relerr_list = [draw_data.mdvX_instant_relerr_list, mdvX_instant_relerr];
        draw_data.mdvS_instant_relerr_list = [draw_data.mdvS_instant_relerr_list, mdvS_instant_relerr];

        iter_data_prev = iter_data;
        iter_data = iter_data_next;
        iter_data_next = three_step_ADMM(iter_data_next, sig_data, SDP_data);
        log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data);
        draw_data = set_draw(iter, draw_data, sig_data, log_data);
        
        if log_data.kkt < tol 
            break;
        end

        if nargin == 6
            sig_data.sig = sig_list(iter);
        end
    end
end

function [x_list, y_list, z_list] = get_Z_3D_lists(mats)
    x_list = zeros(1, length(mats.Z));
    y_list = zeros(1, length(mats.Z));
    z_list = zeros(1, length(mats.Z));
    for i = 1: length(mats.Z)
        Z_mat = mats.Z{i};
        x_list(i) = Z_mat(1, 1);
        y_list(i) = Z_mat(1, 2);
        z_list(i) = Z_mat(2, 2);
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

    % third plot: mdvZ_instant_relerr compare
    rows = 450; cols = 450; linewidth = 1.5;
    % Prepare figure
    f3 = figure('Position', [0, 0, cols, rows], 'Visible', 'on'); set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    % plot
    hold on; grid on; 
    for ii = 1: length(data)
        legend_name = sprintf("$\\log_{10}(t) = %.1f$", log10(data{ii}.t)); 
        plot(log10(min(data{ii}.mdvZ_instant_relerr_list, 1)), 'LineWidth', linewidth, 'Color', LineColors(ii, :), 'DisplayName', legend_name);
    end
    xlabel('Iteration', 'FontSize', fontsize, 'Interpreter', 'latex');
    ylabel('$\log_{10}(\parallel \frac{1}{2} \phi(\bar{Z}; Z^{(k)} - \bar{Z}) - \Delta Z^{(k)} \parallel_F / \parallel \Delta Z^{(k)} \parallel_F)$', 'FontSize', fontsize, 'Interpreter', 'latex');
    legend('show', 'Interpreter', 'latex', 'NumColumns', ceil(length(data)/3), 'Location', 'northoutside', 'FontSize', fontsize);
    % Axes formatting
    ax = gca; set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02)); ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls; ax.GridLineWidth  = 0.6; ax.GridColor = [0, 0, 0]; ax.GridAlpha = 0.4;
    axis tight;
    % Save & close
    print(f3, img_prefix + "taskI_mdvZ_instant_relerr_comp.png", '-dpng', ['-r' num2str(500)]); close(f3);

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
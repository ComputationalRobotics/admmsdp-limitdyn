clc; clear; close all;
cd(fileparts(mfilename('fullpath')));

addpath(genpath("../utils/"));
addpath(genpath("../admm-single-light/"));

%% type
% supported type: 1dc1024, 1et2048, 1tc2048, cancer100, cnhil10, foot, G40mb, G48mb, G48mc, 
% hand, neosfbr25, neosfbr30e8, neu1g, neu2g, neu3g, r12000, swissroll, texture
type_list = [
    "1dc1024", "1et2048", "1tc2048", "cancer100", "cnhil10", "foot", "G40mb", "G48mb", "G48mc", ...
    "hand", "neosfbr25", "neosfbr30e8", "neu1g", "neu2g", "neu3g", "r12000", "swissroll", "texture"
];

for type = type_list

    %% input and output
    fasrc_path = "/n/netscratch/hankyang_lab/Lab/ksc19/2025-field/codes/";
    if exist(fasrc_path, 'dir')
        input_path = fasrc_path + "sig-wave/data/";
    else
        input_path = "./data/";
    end
    input_data_path = input_path + type + "/sig-wave/";
    output_figs_path = "./figs/" + type + "/";
    if ~exist(output_figs_path, "dir")
        mkdir(output_figs_path);
    else
        empty_folder(output_figs_path);
    end

    fprintf("plot sig wave... type: %s \n\n", type);

    %% plot and save
    % options
    options.rows = 300;
    options.cols = 300;
    options.linewidth = 1.2;
    options.ms = 12;
    options.xls = 13;
    options.yls = 13;
    options.fontsize = 13;
    LineColors = linspecer(4, 'qualitative');
    LineColors = flipud(LineColors);
    % load files
    files = dir(fullfile(input_data_path, 'draw_data_Ampl*_N*.mat'));
    for k = 1:numel(files)
        fname = files(k).name;                          % just the file name
        fullpath = fullfile(input_data_path, fname);        % folder + file
        % Extract Ampl and N from the *file name*
        tokens = regexp(fname, 'Ampl=(\d+)_N=(\d+)\.mat', 'tokens', 'once');
        Ampl = str2double(tokens{1});
        N    = str2double(tokens{2});

        % only plot Ampl == 10
        if Ampl == 10
            % Load the file using the full path
            load(fullpath);
            fprintf('Loaded %s: Ampl=%g, N=%d\n', fullpath, Ampl, N);
            img_name = output_figs_path + "four_items.png";
            draw_four_items(draw_data, img_name);
        end

        % % plot diff_norms
        % img_name_1 = output_figs_path + sprintf("diff_norm_Ampl=%d_N=%d.png", ceil(Ampl), N);
        % options.LineColors = LineColors(1:2, :);
        % draw_diff_norm(draw_data, img_name_1, options);
        % % plot anchor_ang
        % img_name_2 = output_figs_path + sprintf("anchor_ang_Ampl=%d_N=%d.png", ceil(Ampl), N);
        % options.LineColors = LineColors(3:4, :);
        % draw_anchor_ang(draw_data, img_name_2, options);
    end
    % plot legends
    options.LineColors = LineColors;
    img_name_3 = output_figs_path + "four_items_legend.png";
    draw_legend(img_name_3, options);

end



%% helper functions
function draw_four_items(draw_data, img_name)
    rows = 300;
    cols = 300;
    linewidth = 1.7;
    ms = 12;
    xls = 13;
    yls = 13;
    fontsize = 13;
    LineColors = linspecer(4, 'qualitative');
    LineColors = flipud(LineColors);

    % Prepare figure
    f = figure('Position', [100, 100, rows, cols], 'Visible', 'on');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    % Data
    sig_list = draw_data.forward.sig_list;
    X_diff_norm_list = draw_data.forward.X_diff_norm_list;
    S_diff_norm_list = draw_data.forward.S_diff_norm_list;
    pinf_list = draw_data.forward.pinf_list;
    dinf_list = draw_data.forward.dinf_list;

    % Log–log plot, two lines with different colors
    plot(log10(sig_list), log10(X_diff_norm_list), 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    plot(log10(sig_list), log10(S_diff_norm_list),  'LineWidth', linewidth, 'Color', LineColors(2, :)); hold on;
    plot(log10(sig_list), log10(pinf_list), 'LineWidth', linewidth, 'Color', LineColors(3, :)); hold on;
    plot(log10(sig_list), log10(dinf_list),  'LineWidth', linewidth, 'Color', LineColors(4, :)); hold on;
    grid on;

    xlabel('$\log_{10}(\sigma)$', 'FontSize', fontsize, 'Interpreter', 'latex');

    % Axes formatting
    ax = gca;
    % set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02));
    % setInsetsXTightYMargin(ax, 0.02);
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.GridLineWidth  = 0.6;
    ax.GridColor      = [0, 0, 0];
    ax.GridAlpha      = 0.4;
    axis tight;

    % Save & close
    print(f, img_name, '-dpng', ['-r' num2str(500)]);
    close(f);
end

function draw_diff_norm(draw_data, img_name, options)
    % Unpack options
    fig_height = options.rows;      % used as pixels for figure size
    fig_width  = options.cols;
    linewidth  = options.linewidth;
    xls        = options.xls;
    yls        = options.yls;
    fontsize   = options.fontsize;
    LineColors = options.LineColors;   % assume at least 2x3

    % Prepare figure
    f = figure('Position', [100, 100, fig_width, fig_height], 'Visible', 'on');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    % Data
    sig_list = [flip(draw_data.backward.sig_list), draw_data.forward.sig_list];
    X_diff_norm_list = [flip(draw_data.backward.X_diff_norm_list), draw_data.forward.X_diff_norm_list];
    S_diff_norm_list = [flip(draw_data.backward.S_diff_norm_list), draw_data.forward.S_diff_norm_list];

    % Log–log plot, two lines with different colors
    loglog(sig_list, X_diff_norm_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    loglog(sig_list, S_diff_norm_list,  'LineWidth', linewidth, 'Color', LineColors(2, :)); hold off;
    grid on;

    % Axes formatting
    ax = gca;
    set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02));
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.GridLineWidth  = 0.6;
    ax.GridColor      = [0, 0, 0];
    ax.GridAlpha      = 0.4;

    % Save & close
    print(f, img_name, '-dpng', ['-r' num2str(500)]);
    close(f);
end

function draw_anchor_ang(draw_data, img_name, options)
    % Unpack options
    fig_height = options.rows;      % used as pixels for figure size
    fig_width  = options.cols;
    linewidth  = options.linewidth;
    xls        = options.xls;
    yls        = options.yls;
    fontsize   = options.fontsize;
    LineColors = options.LineColors;   % assume at least 2x3

    % Prepare figure
    f = figure('Position', [100, 100, fig_width, fig_height], 'Visible', 'on');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    % Data
    sig_list = [flip(draw_data.backward.sig_list), draw_data.forward.sig_list];
    sig_list = sig_list(2:end-1);
    X_anchor_ang_list = [flip(draw_data.backward.X_anchor_diff_ang_list), draw_data.forward.X_anchor_diff_ang_list];
    S_anchor_ang_list = [flip(draw_data.backward.S_anchor_diff_ang_list), draw_data.forward.S_anchor_diff_ang_list];

    % Log–log plot, two lines with different colors
    loglog(sig_list, X_anchor_ang_list, 'LineWidth', linewidth, 'Color', LineColors(1, :)); hold on;
    loglog(sig_list, S_anchor_ang_list,  'LineWidth', linewidth, 'Color', LineColors(2, :)); hold off;
    grid on;

    xlabel('$\sigma$', 'FontSize', fontsize, 'Interpreter', 'latex');

    % Axes formatting
    ax = gca;
    set(ax, 'LooseInset', max(get(ax, 'TightInset'), 0.02));
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.GridLineWidth  = 0.6;
    ax.GridColor      = [0, 0, 0];
    ax.GridAlpha      = 0.4;

    % Save & close
    print(f, img_name, '-dpng', ['-r' num2str(500)]);
    close(f);
end

function draw_legend(img_name, options)
    rows = 50;
    cols = 650;
    linewidth = 3;
    ms = options.ms;
    xls = options.xls;
    yls = options.yls;
    fontsize = options.fontsize;
    LineColors = options.LineColors;

    figure('Position', [0, 0, cols, rows], 'Visible', 'on');
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    legend_1 = "$\log_{10}(\parallel \Delta X^{(k)} \parallel_F)$"; plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    legend_2 = "$\log_{10}(\parallel \Delta S^{(k)} \parallel_F)$"; plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    legend_3 = "$\log_{10}(r_p^{(k)})$"; plot(nan, 'Color', LineColors(3, :), 'LineWidth', linewidth); 
    legend_4 = "$\log_{10}(r_d^{(k)})$"; plot(nan, 'Color', LineColors(4, :), 'LineWidth', linewidth); 

    hold on;
    plot(nan, 'Color', LineColors(1, :), 'LineWidth', linewidth);
    plot(nan, 'Color', LineColors(2, :), 'LineWidth', linewidth);
    plot(nan, 'Color', LineColors(3, :), 'LineWidth', linewidth);
    plot(nan, 'Color', LineColors(4, :), 'LineWidth', linewidth);
    hold off;

    legend({legend_1, legend_2, legend_3, legend_4}, 'Interpreter', 'latex', 'NumColumns', 4, 'Location', 'layout', 'FontSize', fontsize);
    set(gca, 'Visible', 'off');
    
    set(gca, 'LooseInset', [0,0,0,0]) % Remove extra margins
    set(gca, 'Box', 'off') % Remove box around the axes
    set(gcf, 'Color', 'w') % Set background color to white
    axis off % Hide axes

    print(img_name, '-dpng', ['-r', num2str(500)]);
    close all;
end

function setInsetsXTightYMargin(ax, yMargin, xMargin)
    % setInsetsXTightYMargin  Make x tight (left/right), keep y with margins (bottom/top).
    %
    %   setInsetsXTightYMargin()                  uses gca, yMargin=0.04, xMargin=0
    %   setInsetsXTightYMargin(ax)                uses yMargin=0.04, xMargin=0
    %   setInsetsXTightYMargin(ax, yMargin)       sets vertical margin
    %   setInsetsXTightYMargin(ax, yMargin, xMargin) sets both (xMargin usually 0)
    %
    % LooseInset is [left bottom right top] in units of the axes Position.

    if nargin < 1 || isempty(ax), ax = gca; end
    if nargin < 2 || isempty(yMargin), yMargin = 0.04; end
    if nargin < 3 || isempty(xMargin), xMargin = 0; end

    ti = get(ax, 'TightInset');      % [L B R T]
    li = ti;                         % start from tight

    % x margins (left/right) usually 0 to be "tight"
    li(1) = max(ti(1), xMargin);     % left
    li(3) = max(ti(3), xMargin);     % right

    % y margins (bottom/top)
    li(2) = max(ti(2), yMargin);     % bottom
    li(4) = max(ti(4), yMargin);     % top

    set(ax, 'LooseInset', li);
end








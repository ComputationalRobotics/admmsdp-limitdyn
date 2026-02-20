function draw_data = set_draw(iter, draw_data, sig_data, log_data, if_first_time)
    if nargin == 5 && if_first_time == true
        draw_data = struct();
        draw_data.Z_diff_norm_list = []; 
        draw_data.Z_diff_ang_list = []; 
        draw_data.Z_rank_list = [];
        draw_data.Z_norm_list = [];
        draw_data.pinf_list = []; 
        draw_data.dinf_list = []; 
        draw_data.relgap_list = []; 
        draw_data.pobj_list = []; 
        draw_data.dobj_list = []; 
        draw_data.sig_list = [];

        draw_data.X_norm_list = [];
        draw_data.S_norm_list = [];
        draw_data.X_diff_norm_list = [];
        draw_data.S_diff_norm_list = [];
        draw_data.X_diff_ang_list = [];
        draw_data.S_diff_ang_list = [];
    end

    draw_data.Z_diff_norm_list = [draw_data.Z_diff_norm_list, log_data.Z_diff_norm];
    draw_data.Z_diff_ang_list = [draw_data.Z_diff_ang_list, log_data.Z_diff_ang];
    draw_data.Z_rank_list = [draw_data.Z_rank_list, log_data.Z_rank];
    draw_data.pinf_list = [draw_data.pinf_list, log_data.pinf];
    draw_data.dinf_list = [draw_data.dinf_list, log_data.dinf];
    draw_data.relgap_list = [draw_data.relgap_list, log_data.relgap];
    draw_data.pobj_list = [draw_data.pobj_list, log_data.pobj];
    draw_data.dobj_list = [draw_data.dobj_list, log_data.dobj];
    draw_data.sig_list = [draw_data.sig_list, sig_data.sig];
    draw_data.Z_norm_list = [draw_data.Z_norm_list, log_data.Z_norm];

    draw_data.X_norm_list = [draw_data.X_norm_list, log_data.X_norm];
    draw_data.S_norm_list = [draw_data.S_norm_list, log_data.S_norm];
    draw_data.X_diff_norm_list = [draw_data.X_diff_norm_list, log_data.X_diff_norm];
    draw_data.S_diff_norm_list = [draw_data.S_diff_norm_list, log_data.S_diff_norm];
    draw_data.X_diff_ang_list = [draw_data.X_diff_ang_list, log_data.X_diff_ang];
    draw_data.S_diff_ang_list = [draw_data.S_diff_ang_list, log_data.S_diff_ang];
end
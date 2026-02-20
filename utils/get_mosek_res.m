function mosek_data = get_mosek_res(SDP_data)
    At = {SDP_data.At};
    b = SDP_data.b;
    C = {smat_single(SDP_data.C)};
    blk = {['s'], [size(C{1}, 1)]};

    mosek_data = sdpt2mosek(blk, At, C, b);
    mosek_param = [];
    % mosek_param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-12;
    % mosek_param.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-12;
    % mosek_param.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-12;
    % mosek_param.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-14;
    % mosek_param.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL= 1;
    [~, res] = mosekopt('minimize info', mosek_data, mosek_param);
    
    [Xopt, yopt, Sopt, obj] = recover_mosek_sol_blk(res, blk);
    mosek_data.mosek_Xopt_mat = Xopt{1};
    mosek_data.mosek_Sopt_mat = Sopt{1};
    mosek_data.mosek_yopt = yopt;
    mosek_data.mosek_obj = obj;

end
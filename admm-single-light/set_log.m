function log_data = set_log(iter_data_prev, iter_data, iter_data_next, SDP_data)
    % unzip iter_data
    Z_prev = iter_data_prev.Z;
    Z = iter_data.Z;
    Z_next = iter_data_next.Z;
    X = iter_data.X;
    y = iter_data.y;
    S = iter_data.S;
    eigvals = iter_data.eigvals;
    
    % unzip SDP_data
    A = SDP_data.A;
    At = SDP_data.At;
    b = SDP_data.b;
    C = SDP_data.C;
    bscale = SDP_data.bscale;
    Cscale = SDP_data.Cscale;
    objscale = SDP_data.objscale;
    normA = SDP_data.normA;
    normborg = SDP_data.normborg;
    normCorg = SDP_data.normCorg;

    % get KKT
    pinf = norm(normA .* (A * X - b)) * bscale / normborg;
    dinf = norm(At * y + S - C) * Cscale / normCorg;
    pobj = full(C' * X * objscale);
    dobj = full(b' * y * objscale);
    relgap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));

    % get Z log
    Z_rank = nnz(eigvals > 1e-10);
    Z_diff = Z - Z_prev; Z_diff_next = Z_next - Z;
    Z_diff_norm = norm(Z_diff_next);
    Z_norm = norm(iter_data.Z);
    Z_diff_ang = get_ang(Z_diff, Z_diff_next);

    % get X and S log
    X_prev = iter_data_prev.X;
    X_next = iter_data_next.X;
    S_prev = iter_data_prev.S;
    S_next = iter_data_next.S;
    X_norm = norm(X);
    S_norm = norm(S);
    X_diff_norm = norm(X_next - X);
    S_diff_norm = norm(S_next - S);
    X_diff_ang = get_ang(X - X_prev, X_next - X);
    S_diff_ang = get_ang(S - S_prev, S_next - S);
    log_data.X_norm = X_norm;
    log_data.S_norm = S_norm;
    log_data.X_diff_norm = X_diff_norm;
    log_data.S_diff_norm = S_diff_norm;
    log_data.X_diff_ang = X_diff_ang;
    log_data.S_diff_ang = S_diff_ang;

    % collect log data
    log_data.pinf = full(pinf);
    log_data.dinf = full(dinf);
    log_data.pobj = full(pobj);
    log_data.dobj = full(dobj);
    log_data.relgap = full(relgap);
    log_data.Z_rank = full(Z_rank);
    log_data.Z_diff_norm = full(Z_diff_norm);
    log_data.Z_diff_ang = full(Z_diff_ang);
    log_data.kkt = max([pinf, dinf, relgap]);
    log_data.Z_norm = Z_norm;
end
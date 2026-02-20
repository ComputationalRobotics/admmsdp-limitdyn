function iter_data_kpo = three_step_ADMM(iter_data_k, sig_data, SDP_data)
    A = SDP_data.A;
    At = SDP_data.At;
    b = SDP_data.b;
    C = SDP_data.C;
    R = SDP_data.R;
    P = SDP_data.P;

    X_k = iter_data_k.X;
    S_k = iter_data_k.S;

    sig = sig_data.sig;

    rhsy = 1/sig * b - A * (1/sig * X_k + S_k - C);
    y_kpo = chol_solve(rhsy, R, P);
    Xb = X_k + sig * (At * y_kpo - C); 
    if isfield(SDP_data, 'if_gpu_proj') && SDP_data.if_gpu_proj == true
        [Xb_proj, eigvals] = proj_single_gpuarray(Xb);
    else
        [Xb_proj, eigvals] = proj_single(Xb);
    end
    S_kpo = 1/sig * (Xb_proj - Xb);
    X_kpo = X_k + sig * (S_kpo + At * y_kpo - C);

    iter_data_kpo.X = X_kpo;
    iter_data_kpo.y = y_kpo;
    iter_data_kpo.S = S_kpo;
    iter_data_kpo.Z = X_kpo - sig * S_kpo;
    iter_data_kpo.eigvals = eigvals;
    iter_data_kpo.cnt = iter_data_k.cnt + 1;
end
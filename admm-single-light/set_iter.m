function iter_data = set_iter(Z, sig_data, SDP_data, cnt)
    sig = sig_data.sig;
    X = proj_single(Z);
    S = (X - Z) / sig;

    iter_data.X = X;
    iter_data.S = S;
    iter_data.Z = Z;
    iter_data.eigvals = eig(smat_single(Z));

    if nargin == 3
        iter_data.y = zeros(size(SDP_data.b));
    end

    if nargin == 4
        iter_data.cnt = cnt;
    else 
        iter_data.cnt = 1;
    end
end
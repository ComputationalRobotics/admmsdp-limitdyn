function [X_svec, eigvals, eigvecs] = proj_single_gpuarray(X_svec)
    A = smat_single(X_svec);

    % 1) Move to GPU and symmetrize (important for numerical stability)
    Ag = gpuArray(A);
    
    % 2) Eigen-decomposition on GPU
    [V,D] = eig(Ag);         % V and D are on the GPU
    
    % 3) Zero out negative eigenvalues (PSD cone projection)
    lambda = diag(D);                % eigenvalues (on GPU)
    lambda_proj = max(real(lambda), 0);  % clip negatives to 0
    
    D_proj = diag(lambda_proj);
    Ag_psd = V * D_proj * V';        % still on GPU
    
    % 4) Pull back to CPU
    A_psd = gather(Ag_psd);

    X_svec = svec_single(A_psd);

    if (nargout == 2)
        eigvals = gather(lambda);
    elseif (nargout == 3)
        eigvals = gather(lambda);
        eigvecs = gather(V);
    end
end
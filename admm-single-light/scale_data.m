function [iter_data, SDP_data] = scale_data(iter_data, sig_data, SDP_data, if_scale)
    % upzip data
    if_y_exist = isfield(iter_data, 'y');
    A = SDP_data.A; 
    At = SDP_data.At; 
    b = SDP_data.b;
    C = SDP_data.C; 
    X = iter_data.X; 
    if if_y_exist
        y = iter_data.y; 
    end
    S = iter_data.S;
    
    % save b and C's norm
    borg = b;
    Corg = C;
    normborg = 1 + norm(borg);
    normCorg = 1 + norm(Corg);
    
    % scale At
    if nargin == 4 && if_scale == false
        normA = 1;
    else
        normA = sum(At .* At, 1)';
        normA = max(1.0, sqrt(normA));
        DA = spdiags(1./normA, 0, size(At, 2), size(At, 2));
        At = At * DA;
        A = At';
        b = b ./ normA;
        if if_y_exist
            y = normA .* y;
        end
    end
    AAt = A * At;
    if nnz(AAt) > 0.1 * size(AAt, 1)^2
        AAt = full(AAt);
        R = chol(AAt);
        P = speye(size(AAt));
    else
        AAt = sparse(AAt);
        [R, ~, P] = chol(AAt);
    end

    % scale b, C
    normb = 1 + norm(b);
    normC = 1 + norm(C);
    if nargin > 3 && if_scale == false
        bscale = 1;
        Cscale = 1;
    else
        bscale = normb;
        Cscale = normC;
    end
    objscale = bscale * Cscale;
    b = b / bscale;
    C = C / Cscale;
    X = X / bscale;
    S = S / Cscale;
    if if_y_exist
        y = y / Cscale; 
    end

    % save iter_data
    iter_data.X = X;
    if if_y_exist
        iter_data.y = y;
    end
    iter_data.S = S;
    iter_data.Z = X - sig_data.sig * S;

    % save SDP_data
    SDP_data.A = A;
    SDP_data.At = At;
    SDP_data.b = b;
    SDP_data.C = C;
    SDP_data.R = R;
    SDP_data.P = P;
    SDP_data.bscale = bscale;
    SDP_data.Cscale = Cscale;
    SDP_data.objscale = objscale;
    SDP_data.normA = normA;
    SDP_data.normborg = normborg;
    SDP_data.normCorg = normCorg;
end




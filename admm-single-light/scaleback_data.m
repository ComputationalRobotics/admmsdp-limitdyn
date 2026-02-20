function rescaled_data = scaleback_data(iter_data, sig_data, SDP_data)
    bscale = SDP_data.bscale;
    Cscale = SDP_data.Cscale;
    normA = SDP_data.normA;
    sig = sig_data.sig;

    X = iter_data.X;
    y = iter_data.y;
    S = iter_data.S;
    rescaled_data.X = bscale * X;
    rescaled_data.y = Cscale * (y ./ normA);
    rescaled_data.S = Cscale * S;
    rescaled_data.Z = rescaled_data.X - sig * rescaled_data.S;
    rescaled_data.cnt = iter_data.cnt;
end
function [iter_data, SDP_data, sig_data] = rescale_Candb(iter_data, sig_data, SDP_data, Cscale_inc, bscale_inc)
    % upzip data
    b = SDP_data.b;
    C = SDP_data.C; 
    bscale = SDP_data.bscale;
    Cscale = SDP_data.Cscale;
    objscale = SDP_data.objscale;
    X = iter_data.X; 
    S = iter_data.S;
    y = iter_data.y;
    sig = sig_data.sig;
    
    % rescaling b and C
    bscale = bscale * bscale_inc;
    Cscale = Cscale * Cscale_inc;
    objscale = objscale * bscale_inc * Cscale_inc;
    b = b ./ bscale_inc;
    C = C ./ Cscale_inc;
    X = X ./ bscale_inc;
    S = S ./ Cscale_inc;
    y = y ./ Cscale_inc;
    sig = sig * Cscale_inc / bscale_inc;
    
    % pass rescaled data back
    sig_data = set_sig(sig);
    SDP_data.b = b;
    SDP_data.C = C;
    SDP_data.bscale = bscale;
    SDP_data.Cscale = Cscale;
    SDP_data.objscale = objscale;
    iter_data.X = X;
    iter_data.y = y;
    iter_data.S = S;
    iter_data.Z = X - sig_data.sig * S;
end
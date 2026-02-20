function sig_data = set_sig(sig_init, sig_data)
    sig_data.sig = sig_init;
    sig_data.prim_win = 0;
    sig_data.dual_win = 0;

    if nargin == 2
        sig_data.stop_iter = sig_data.stop_iter;
    end
end
function sig_data = update_sig_2(iter, sig_data, log_data)
    sig = sig_data.sig;
    prim_win = sig_data.prim_win;
    dual_win = sig_data.dual_win;
    
    pinf = log_data.pinf;
    dinf = log_data.dinf;

    ratioconst = 1e0;
    feasratio = ratioconst * pinf / dinf;

    if (feasratio < 1)
        prim_win = prim_win + 1;
    else
        dual_win = dual_win + 1;
    end

    sigmax = 1e6;
    sigmin = 1e-6;
    sigscale = max(1.0, 1.5 * exp(-iter/4000));

    if ( ( iter > 100 && iter <= 500 && rem(iter, 100) == 1 ) ...
          || ( iter > 500 && iter <= sig_data.stop_iter && rem(iter, 200) == 1 )  )
        if (prim_win > 1.35 * dual_win)
            prim_win = 0;
            sig = min(sigmax, sig * sigscale);
        elseif (dual_win > 1.35 * prim_win)
            dual_win = 0;
            sig = max(sigmin, sig / sigscale);
        end
    end

    sig_data.sig = sig;
    sig_data.prim_win = prim_win;
    sig_data.dual_win = dual_win;
end

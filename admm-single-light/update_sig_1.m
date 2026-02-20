function sig_data = update_sig_1(iter, sig_data, log_data)
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

    sigmax = 1e3;
    sigmin = 1e-3;
    % sigscale = max(1.01, 2 * exp(-iter/50000));
    sigscale = 1.2;
    
    if ~isfield(sig_data, 'stop_iter')
        sig_data.stop_iter = 2e4;
    end
    if iter > sig_data.stop_iter 
        return;
    end

    if ( (iter <= 200 && rem(iter, 50) == 1) ...
        || (iter > 200 && iter <= 1000 && rem(iter, 100) == 1) ...
        || (iter > 1000 && iter <= 5000 && rem(iter, 200) == 1) ...
        || (iter > 5000 && rem(iter, 1000) == 1) )
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

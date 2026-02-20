function sig_data = update_sig_wave(iter, sig_data, log_data)
    if ~isfield(sig_data, 'a')
        a = 0;
    else 
        a = sig_data.a;
    end

    if ~isfield(sig_data, 'b')
        b = 6/24;
    else
        b = sig_data.b;
    end
    
    if ~isfield(sig_data, 'N')
        N = ceil(1e5);
    else
        N = sig_data.N;
    end
    
    if ~isfield(sig_data, 'Ampl')
        Ampl = 1e4;
    else
        Ampl = sig_data.Ampl;
    end

    if ~isfield(sig_data, 'sig_base')
        sig_base = 1.0;
    else
        sig_base = sig_data.sig_base;
    end

    x = mod(iter, N) / N;
    res = wave_1(x, a, b);
    sig_data.sig = sig_base * Ampl^res;
    
    sig_data.a = a;
    sig_data.b = b;
    sig_data.N = N;
    sig_data.Ampl = Ampl;
    sig_data.sig_base = sig_base;
    sig_data.wave_func = @(x) wave_1(x, a, b);
end


function res = wave_1(x, a, b)
    assert(x >= 0 && x <= 1);
    assert(a >= 0 && b >= 0 && a+b <= 0.25);

    function res = f(x, a, b)
        % x should be [0, 0.25]!
        if x >= 0 && x <= a 
            res = 0;
        elseif x > a && x < a+b 
            res = 1/b * x - a/b;
        elseif x >= a+b && x <= 0.25 
            res = 1;
        else
            error("x should be in [0, 0.25]! \n");
        end
    end

    if x <= 0.25 
        res = f(x, a, b);
    elseif x > 0.25 && x <= 0.5 
        res = f(0.5-x, a, b);
    elseif x > 0.5 && x <= 0.75 
        res = -f(x-0.5, a, b);
    else
        res = -f(1-x, a, b);
    end
end
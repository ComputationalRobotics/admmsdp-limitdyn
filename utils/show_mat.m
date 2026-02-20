function f = show_mat(A_mat, if_log)
    Amp = A_mat ./ max(abs(A_mat(:)));
    if nargin == 2 && if_log == true
        Amp = abs(Amp);
        eps = 1e-16;
        Amp = log10( Amp./max(Amp(:)) + eps );   % normalize & avoid -Inf
    end
    
    f = figure;
    if nargin == 2 && if_log == true
        imagesc(Amp, [-8 0]);        % show -80 to 0 dB
    else
        imagesc(Amp);
    end
    axis image;
    colormap(parula);
    colorbar;
end
function surf_mat(A_mat, if_log)
    Amp = A_mat ./ max(abs(A_mat(:)));
    if nargin == 2 && if_log == true
        Amp = abs(Amp);
        eps = 1e-16;
        Amp = log10( Amp./max(Amp(:)) + eps );   % normalize & avoid -Inf
    end

    figure;
    surf(Amp,'EdgeColor','none');   % no mesh lines
    shading interp;
    colormap(turbo);
    colorbar;
    view(45,40);
end
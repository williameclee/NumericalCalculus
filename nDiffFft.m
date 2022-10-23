function F_d1 = nDiffFft(F,varargin)
%% nDiffFft Differentiation with fast Fourier transformation 
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - order: The order of derivative. The default value is 1 (first 
%     derivative).
%   - overlap: The amount of overlapping grids if the arrray is periodic.
%     The default value is 0 (no overlap).
% OUTPUT:
%   - F_d1: The first derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFft(Y,Dx,"order",1,"overlap",1)

    h       = 1;
    order   = 1;
    overlap = 0;

    p = inputParser;
    addRequired(p,"F");
    addOptional(p,"h",h);
    addOptional(p,"order",order);
    addOptional(p,"overlap",overlap);

    F = reshape(F,[1,length(F)]);

    N = length(F)-overlap;
    F = F(1:N);
    L = h*(N+1);
    
    if mod(length(F),2) == 0
        omega = fftshift(1i*(2*pi/L)*(-(N/2):(N/2)-1));
    else
        omega = fftshift(1i*(2*pi/L)*(-floor(N/2):floor(N/2)));
    end
    F_hat    = fft(F);
    F_d1_hat = omega.^order.*F_hat;

    if isreal(F)
        F_d1 = real(ifft(F_d1_hat));
    end
    F_d1 = [F_d1,F_d1(1:overlap)];
end
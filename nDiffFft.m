function Fd = nDiffFft(F,varargin)
%% nDiffFft Differentiation with fast Fourier transformation 
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - dimension: The dimsion along which the derivative of F is calculated. 
%   - order: The order of derivative. The default value is 1 (first 
%     derivative).
%   - overlap: The amount of overlapping grids if the arrray is periodic.
%     The default value is 0 (no overlap).
%   - boundary: The spectral method in its natural is periodic. This 
%     argument is only included for continuity.
% OUTPUT:
%   - Fd: The first derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd4(Y,Dx,"dimension",2,"order",1, ...
%                    "overlap",0,"boundary","periodic")

    % Assigning default parameters
    h        = 1;
    dim      = 1;
    order    = 1;
    overlap  = 0;
    boundary = "periodic";
    % Updating parameters
    p = inputParser;
        addRequired(p,"F");
        addOptional(p,"h",h);
        addOptional(p,"dimension",dim);
        addOptional(p,"order",order);
        addParameter(p,"overlap",overlap);
        addParameter(p,"boundary",boundary);
        parse(p,F,varargin{:});
    F        = p.Results.F;
    h        = p.Results.h;
    dim      = p.Results.dimension;
    order    = p.Results.order;
    overlap  = p.Results.overlap;
    % Getting rid of overlapping grids
    id      = repmat({':'},1,ndims(F));
    id(dim) = {1:size(F,dim)-overlap};
    F       = F(id{:});
    % Reshaping the array if necessary
    F_size = size(F);
    if ismatrix(F) && F_size(1) == 1
        F = reshape(F,[length(F),1]);
    end
    N = size(F,dim);
    L = h*(N+1);
    
    if mod(N,2) == 0
        omega = fftshift(1i*(2*pi/L)*(-(N/2):(N/2)-1));
    else
        omega = fftshift(1i*(2*pi/L)*(-floor(N/2):floor(N/2)));
    end
    id      = F_size;
    id(dim) = 1;
    omega = repmat(omega,id);
    F_hat  = fft(F,[],dim);
    Fd_hat = omega.^order.*F_hat;
    Fd     = ifft(Fd_hat,[],dim);
    % Getting rid of imaginary parts if necessary
    if isreal(F)
        Fd = real(Fd);
    end
    Fd = reshape(Fd,F_size);
    % Adding overlapping grids back
    id      = repmat({':'},1,ndims(F));
    id(dim) = {1:overlap};
    Fd = cat(dim,Fd,Fd(id{:}));
end
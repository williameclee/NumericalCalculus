function Fd = nDiffFd4(F,varargin)
%% nDiffFd4 Differentiation with four-point (central) finite difference
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - order: The order of derivative. The default value is 1 (first 
%     derivative); highest supported value is 2.
%   - "boundary": The boundary conditions at ends of the array. If it is
%     not "periodic", forward and backward differences will be utilised to
%     achieve the same order of truincation error. The default value is
%     "wall".
% OUTPUT:
%   - Fd: The derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd4(Y,Dx,"order",2,"boundary","periodic")

    h     = 1;
    order = 1;
    boundary = "wall";
    
    p = inputParser;
        addRequired(p,'F');
        addOptional(p,'h',h);
        addOptional(p,'order',order);
        addParameter(p,'boundary',boundary);
        parse(p,F,varargin{:});
    F        = p.Results.F;
    h        = p.Results.h;
    order    = p.Results.order;
    boundary = p.Results.boundary;

    F = reshape(F,[1,length(F)]);

    if order == 1 % first derivative
        Fd = (+ 8 *circshift(F,-1) - 8*circshift(F,1) ...
              - 1 *circshift(F,-2) + 1*circshift(F,2))/(12*h);
        if ~strcmp(boundary,"periodic")
            Fd(1:2)       = (- 25*F(1:2) + 48*F(2:3) ...
                             - 36*F(3:4) + 16*F(4:5) ...
                             -  3*F(5:6))/(12*h);
            Fd(end-1:end) = (+ 25*F(end-1:end)   - 48*F(end-2:end-1) ...
                             + 36*F(end-3:end-2) - 16*F(end-4:end-3) ...
                             +  3*F(end-5:end-4))/(12*h);
        end
    elseif order == 2 % secend derivative
        Fd = (- 30*F ...
              + 16*circshift(F,-1) + 16*circshift(F,1) ...
              -  1*circshift(F,-2) -  1*circshift(F,2))/(12*h^2);
        if ~strcmp(boundary,"periodic")
            Fd(1:2)       = (+  225*F(1:2) - 770*F(2:3) ...
                             + 1070*F(3:4) - 780*F(4:5) ...
                             +  305*F(5:6) -  50*F(6:7))/(60*h^2);
            Fd(end-1:end) = (+  225*F(end-1:end)   - 770*F(end-2:end-1) ...
                             + 1070*F(end-3:end-2) - 780*F(end-4:end-3) ...
                             +  305*F(end-5:end-4) -  50*F(end-6:end-5))/(60*h^2);
        end
    end
end
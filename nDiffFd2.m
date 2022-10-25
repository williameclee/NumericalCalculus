function Fd = nDiffFd2(F,varargin)
%% nDiffFd2 Differentiation with two-point (central) finite difference
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - dimension: The dimsion along which the derivative of F is calculated. 
%   - order: The order of derivative. The default value is 1 (first 
%     derivative); highest supported value is 2.
%   - overlap: The amount of overlapping grids if the arrray is periodic.
%     The default value is 0 (no overlap).
%   - boundary: The boundary conditions at ends of the array. If it is not 
%     "periodic", forward and backward differences will be utilised to
%     achieve the same order of truincation error. The default value is
%     "wall".
% OUTPUT:
%   - Fd: The derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd2(Y,Dx,"dimension",2,"order",2, ...
%                    "overlap",0,"boundary","periodic")
    
    % Assigning default parameters
    h        = 1;
    dim      = 1;
    order    = 1;
    overlap  = 0;
    boundary = "wall";
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
    boundary = p.Results.boundary;
    % Getting rid of overlapping grids
    id      = repmat({':'},1,ndims(F));
    id(dim) = {1:size(F,dim)-overlap};
    F       = F(id{:});
    % Reshaping the array if necessary
    F_size = size(F);
    if ismatrix(F) && F_size(1) == 1
        F = reshape(F,[length(F),1]);
    end
    % Calculating the derivatives
    if order == 1 % first derivative
        % Central
        Fd = (+ 1*circshift(F,-1,dim) - 1*circshift(F,+1,dim))/(2*h);
        if ~(strcmp(boundary,"periodic") || strcmp(boundary,"p"))
            % Forward
            id      = repmat({':'},1,ndims(F));
            id(dim) = {1:3};
            F_f  = F(id{:});
            Fd_f = (- 3*circshift(F_f,-0,dim) + 4*circshift(F_f,-1,dim) ...
                    - 1*circshift(F_f,-2,dim))/(2*h);
            id(dim) = {1};
            Fd(id{:}) = Fd_f(id{:});
            % Backward
            id      = repmat({':'},1,ndims(F));
            id(dim) = {size(F,dim)-2:size(F,dim)};
            F_e  = F(id{:});
            Fd_e = (+ 3*circshift(F_e,+0,dim) - 4*circshift(F_e,+1,dim) ...
                    + 1*circshift(F_e,+2,dim))/(2*h);
            id(dim) = {size(F,dim)};
            id_e    = id; id_e(dim) = {size(Fd_e,dim)};
            Fd(id{:}) = Fd_e(id_e{:});
        end
    elseif order == 2 % secend derivative
        % Central
        Fd = (+ 1*circshift(F,-1,dim) - 2*F + 1*circshift(F,+1,dim))/(1*h^2);
        if ~(strcmp(boundary,"periodic") || strcmp(boundary,"p"))
            % Forward
            id      = repmat({':'},1,ndims(F));
            id(dim) = {1:3};
            F_f  = F(id{:});
            Fd_f = (+ 1*circshift(F_f,-0,dim) - 2*circshift(F_f,-1,dim) ...
                    + 1*circshift(F_f,-2,dim))/(1*h^2);
            id(dim) = {1};
            Fd(id{:}) = Fd_f(id{:});
            % Backward
            id      = repmat({':'},1,ndims(F));
            id(dim) = {size(F,dim)-2:size(F,dim)};
            F_e  = F(id{:});
            Fd_e = (+ 1*circshift(F_e,+0,dim) - 2*circshift(F_e,+1,dim) ...
                    + 1*circshift(F_e,+2,dim))/(1*h^2);
            id(dim) = {size(F,dim)};
            id_e    = id; id_e(dim) = {size(Fd_e,dim)};
            Fd(id{:}) = Fd_e(id_e{:});
        end
    end
    Fd = reshape(Fd,F_size);
    % Adding overlapping grids back
    id      = repmat({':'},1,ndims(F));
    id(dim) = {1:overlap};
    Fd = cat(dim,Fd,Fd(id{:}));
end
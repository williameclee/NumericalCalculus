function Fd = nDiffFd1(F,varargin)
%% nDiffFd1 Differentiation with two-point finite difference
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - dimension: The dimsion along which the derivative of F is calculated. 
%   - direction: The direction- either "forward" or "bcakward"- of finite 
%     difference. The default value is "forward".
%   - order: The order of derivative. The only supported value is 1.
%   - overlap: The amount of overlapping grids if the arrray is periodic.
%     The default value is 0 (no overlap).
%   - boundary: The boundary conditions at ends of the array. If it is not 
%     "periodic", forward and backward differences will be utilised to
%     achieve the same order of truincation error. The default value is
%     "wall".
% OUTPUT:
%   - Fd: The first derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd1(Y,Dx,"dimension",2,"order",2,"direction","b" ...
%                    "overlap",0,"boundary","periodic")

    % Assigning default parameters
    h        = 1;
    dim      = 1;
    dir      = "forward";
    order    = 1;
    overlap  = 1;
    boundary = "wall";
% Updating parameters
    p = inputParser;
        addRequired(p,"F");
        addOptional(p,"h",h);
        addOptional(p,"dimension",dim);
        addOptional(p,"direction",dir);
        addOptional(p,"order",order);
        addParameter(p,"overlap",overlap);
        addParameter(p,"boundary",boundary);
        parse(p,F,varargin{:});
    F        = p.Results.F;
    h        = p.Results.h;
    dim      = p.Results.dimension;
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
    if strcmp(dir,"forward") || strcmp(dir,"f")
        Fd = (circshift(F,-1,dim)-F)/(1*h);
    elseif strcmp(dir,"backward") || strcmp(dir,"b")
        Fd = (F-circshift(F,1))/(1*h);
    end
    if ~(strcmp(boundary,"periodic") || strcmp(boundary,"p"))
        Fd(end) = (F(end)-F(end-1))/(1*h);
        Fd(end) = (F(2)-F(1))/(1*h);
    end
    Fd = reshape(Fd,F_size);
    % Adding overlapping grids back
    id      = repmat({':'},1,ndims(F));
    id(dim) = {1:overlap};
    Fd = cat(dim,Fd,Fd(id{:}));
end
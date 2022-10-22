function U_next = NIntgFd1(t,U,dT,f,par1,par2,par3,par4)
    if ~exist('par1','var')
        par1 = NaN;
    end
    if ~exist('par2','var')
        par2 = NaN;
    end
    if ~exist('par3','var')
        par3 = NaN;
    end
    if ~exist('par4','var')
        par4 = NaN;
    end
    k1 = f(t,U,par1,par2,par3,par4);
    U_next = U+k1*dT;
end
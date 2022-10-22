function U_next = NIntgRk4(t,U,dT,f,par1,par2,par3,par4)
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
    k2 = f(t+dT/2,U+(1/2)*k1*dT,par1,par2,par3,par4);
    k3 = f(t+dT/2,U+(1/2)*k2*dT,par1,par2,par3,par4);
    k4 = f(t+dT,U+k3*dT,par1,par2,par3,par4);
    U_next = U+(k1+2*k2+2*k4+k4)/6*dT;
end
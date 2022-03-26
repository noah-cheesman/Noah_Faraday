function [omega,omega_s]=omegas(pars,i)
    if length(pars.T)>1
        omega_s=pars.omega(i)*pars.T(i)/(pars.T(i)+pars.c_theta(i));
    else
        omega_s=pars.omega*pars.T/(pars.T+pars.c_theta);
    end
	if length(pars.T)>1
        omega=pars.omega(i);
    else
        omega=pars.omega;
    end
end
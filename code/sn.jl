include("nu.jl")

type sn
  max_date::Float64
  ra::Float64
  dec::Float64
  z::Float64
  zenith_bin::Int
  nb::Float64
  coefs::Array{Float64,1}

  function sn(max_date::Real,ra::Float64,dec::Float64,z::Float64)
      new_sn = new(max_date,ra,dec,z);
      new_sn.coefs = Float64[];

      return new_sn
  end
end

add_coefs!(t_sn::sn, coefs::Array{Float64,1}) = t_sn.coefs=coefs;
add_nb_and_zenith_bin_idx!(t_sn::sn, nb::Float64, kk::Int) = (t_sn.nb=nb; t_sn.zenith_bin = kk;)

rm_coefs!(t_sn) = t_sn.coefs=Float64[];

# Function that determines which nus from an input array are associated with a
# single sn. This is done by checking to see if the sn is within the error
# region defined by its uncertainty in arrival direction and time
#
# Parameters:
#   alpha   -- the size of uncertainty window in arrival direction  (0,1]
#   n_hi    -- the largest amount of days before sn Max Date that neutrinos will
#              will be accepted in the window
#   n_lo    -- the smallest amount of days before sn Max Date that neutrinos will
#              will be accepted in the window
#   t_coef  -- coefficient used in the Taylor series approximation of the Kent
#              Distribution
#
# Input :
#   t_sn  -- a single sn object
#   t_nus -- an array of nu objects
#
# Output :
#   an array of Bools indicating which of the neutrinos are associated with
#   the input sn object

function find_associated_nus(t_sn::sn, t_nus::Array{nu,1})
    n_hi = 19;
    n_lo = 4;
    alpha = 0.99;
    t_coef = 1.0-2.0*alpha;

    t_len_nu = length(t_nus);

    # determine if the neutrino is within the time window of the potential
    # core collapse

    in_time_window = [t_sn.max_date - n_hi <= t_nus[j].mjd <= t_sn.max_date - n_lo for j in 1:t_len_nu];

    # Parameters used in the Kent Distribution

    t_kappa = [1./t_nus[j].ang_err.^2 for j in 1:t_len_nu];
    t_mu = [sin(t_sn.dec)*sin(t_nus[j].dec) + cos(t_sn.dec)*
        cos(t_nus[j].dec)*cos(t_sn.ra-t_nus[j].ra) for j in 1:t_len_nu];

    in_ang_window = Array(Bool,t_len_nu);
    acceptance_mu = Array(Float64,t_len_nu);

    # determine if the neutrino is within the angular acceptance window of
    # the sn. Here we utilize different approximations of the
    # Kent Distribution 

    for j in 1:t_len_nu
        if t_kappa[j] > 10.0
            acceptance_mu[j] = 1.0 + log(1.0-alpha)/t_kappa[j];
        elseif 1e-2 < t_kappa[j] <= 10.0
            acceptance_mu[j] = log(exp(t_kappa[j]) - 2.0*alpha*sinh(t_kappa[j]))/t_kappa[j];
        elseif 0.0 < t_kappa[j] <= 1e-2
            acceptance_mu[j] = log(1+t_coef*t_kappa[j] + 0.5*t_kappa[j]^2 + t_coef*t_kappa[j]^3/6)/t_kappa[j]
        else
            error("strange kappa in find_associated_nus");
        end
    end
    in_ang_window[:] = t_mu .> acceptance_mu;
    return in_time_window.*in_ang_window;
end

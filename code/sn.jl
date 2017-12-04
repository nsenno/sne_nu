include("nu.jl")

type sn
  max_date::Float64
  ra::Float64
  dec::Float64
  z::Float64
  zenith_bin::Int
  nb::Float64
  associated_nus::Array{nu,1}
  coefs::Array{Float64,1}

  function sn(max_date::Real,ra::Float64,dec::Float64,z::Float64)
      new_sn = new(max_date,ra,dec,z);
      new_sn.associated_nus=nu[];
      new_sn.coefs = Float64[];

      return new_sn
  end
end

add_coefs!(t_sn::sn, coefs::Array{Float64,1}) = t_sn.coefs=coefs;
add_nb!(t_sn::sn, nb::Float64) = t_sn.nb=nb;
add_nu!(t_sn::sn, nu::Array{nu,1}) = append!(t_sn.associated_nus,nu);
add_nu!(t_sn::sn, nu::nu) = append!(t_sn.associated_nus,nu);

rm_associated_nus!(t_sn) = t_sn.associated_nus=nu[];
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

# function that calculates the average number of neutrinos associated with each SN using scrambled neutrino data
# the results are written to a file "nb_data.dat"
function calc_nb(N_samples = 1000)

    len_sne = length(sne);
    len_nu = length(nu);

    # Array to store the number of background neutrinos for each SN for each randomized sample

    num_nb = zeros(Int,N_samples,len_sne);

    for n in 1:N_samples
        sample_nus = create_nu_array(create_sample_nus(nu_data,len_nu))
        num_nb[n,:] = [sum(find_associated_nus(sn,sample_nus)) for sn in sne]
    end

    writedlm(string(data_dir,"nb_data.dat"),[mean(num_nb[:,i]) for (i,_) in enumerate(sne)])

end

# Function that calculates the coefficients neccesary to simplify the TS
# calculation. In this version, only arrival time and direction are considered.
#
# Input :
#   t_sn  --  a single sn object that already has its associated nus
#             array filled with nu objects
#
# Output :
#   None 

function calc_coefs!(t_sn::sn)
    if length(t_sn.associated_nus) > 0

        # Signal direction PDF from Kent distirbution

        temp_s = map(x-> S_dir(t_sn,x),t_sn.associated_nus);

        # Signal time PDF from Poisson distribution with mean 13 days

        temp_s .*= map(x-> S_time(t_sn,x),t_sn.associated_nus);

        # Background direction PDF comprized of RA (uniform in [0,2pi]) and
        # Dec. which is computed experimentally from the function "get_dec_pdf"
        # from "sig_and_bkg_pdfs.jl"

        temp_b = (1/(2*pi))*map(B_dec,[t_sn.associated_nus[j].dec for j in 1:length(t_sn.associated_nus)]);

        # Background time PDF is uniform over 16 days (which is the size
        # of the acceptance window)

        temp_b ./= 16.0;

        # Add the relevant coeficient to input sn

        add_coefs!(t_sn,temp_s./temp_b./t_sn.nb);
    end
end

using KernelDensity

include("sn.jl")

# Function that calculates the signal probability distribution for an
# individual neutrino with respect to an individual SN. We use the Kent
# Distribution to account for the spherical coordinates. The angular uncertainty
# of the neutrino is also taken into account. Because of numerical
# effects, the Taylor series approximation is used for very large and very
# small values of kappa (i.e., the angular uncertainty)
#
# Inputs :
#    t_sn  --  an individual sn object to be compared with t_nu
#    t_nu  --  an individual nu object to be compared with t_sn
#
# Outputs :
#    the Kent Distribution PDF evaluated using the angular distance between the
#    neutrino and sn, accounting for the neutrino angular uncertainty
#
# Parameters :
#    mu    --  the cosine of the angular distance between the t_sn and t_nu
#    kappa --  the inverse of the angular variance (only due to the neutrino)

function S_dir(t_sn::sn,t_nu::nu)

    mu = sin(t_sn.dec)*sin(t_nu.dec) + cos(t_sn.dec)*cos(t_nu.dec)*cos(t_sn.ra-t_nu.ra);
    kappa = 1/t_nu.ang_err^2;

    if kappa > 10.0

        result = kappa/(2*pi)*exp(kappa*(mu-1.0))

    elseif 1e-2 < kappa <= 10.0

        result = kappa/(4*pi*sinh(kappa))*exp(kappa*mu);

    elseif 0.0 < kappa <= 1e-2

        result = (1.0 + kappa*mu)/(4*pi);

    else

        error("strange kappa");

    end

    @assert result >= 0.0;

    return cos(t_nu.dec)*result;

end

# Function that calculates the signal PDF for an individual neutrino with
# respect to an individual SN. We use the Poisson distribution assuming an
# mean of 13 days between the maximum SN brightness and the SN core-collapse
# (i.e., when the neutrinos are produced).
#
# Inputs :
#    t_sn  --  an individual sn object to be compared with t_nu
#    t_nu  --  an individual nu object to be compared with t_sn
#
# Outputs :
#    the Poisson PDF evaluated using the difference in time  (in days) between
#    the neutrino detection and sn maximum optical brightness
#

function S_time(t_sn::sn, t_nu::nu)

    # 0.5 accounts for differences in convention between the SNe dataset
    # and the neutrino dataset

    return pdf(Poisson(13.0),t_sn.max_date-t_nu.mjd-0.5)

end

# Function that converts the Dec. values from background neutrinos into an
# an anonymous function that caluculates the background PDF for Dec.
#
# Input :
#    t_nu_dec  --  an array of Dec's from background neutrinos
#
# Output :
#    Anonymous Function :

#       Input :
#          xx   --  Neutrino declination
#       Output :
#          B_dec(xx) if result is positive, or NaN otherwise
#

function get_dec_pdf(t_nu_dec::Array{Float64,1})

    nu_dec_kde = kde(t_nu_dec);
    nu_dec_interp = InterpKDE(nu_dec_kde);

    function (xx) yy = pdf(nu_dec_interp,xx);  return yy > 0 ? yy : NaN end
end

# Read in values of Dec for background neutrinos. Relevant file is created by
# "generate_nu_dec_bkg_pdf.jl"

bkg_dec = readdlm(string(data_dir,"bkg_nu_dec_values.dat"))[:];

# Convert background neutrino Dec values to Background PDF

B_dec = get_dec_pdf(bkg_dec)

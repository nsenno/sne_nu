# Create an immutable type nu which contains the arrival time and direction (Right Assension RA and Declination Dec)
# These objects will be useful when creating random samples

immutable nu
    mjd::Float64
    ang_err::Float64
    ra::Float64
    dec::Float64
end

# Function that reads in an array of neutrino data and converts it into an array of nu objects
#
# Inputs :
#    nu_data --  a 2D array of floats containing raw neutrino data
#                data organization : [Arrival Day (MJD), Angular Error, RA, Dec]
#
# Output :
#    an array of nu objects created from raw neutrino data
#

function create_nu_array(nu_data::Array{Float64,2})
    return [nu(nu_data[i,:]...) for i in 1:length(nu_data[:,1])]
end

# Function that reads in raw neutrino data an returns a randomized sample
#
# Inputs :
#   orig_data --  a 2D array of floats containing original neutrino data
#                 data organization : [Arrival Day (MJD), Angular Error, RA, Dec]
#
#   N_nus     --  Number of randomized samples. If no value is given use the entire length of the orig_data
#
# Output :
#    a 2D array of floats containing a randomized sample of neutrino data
#    data organization : [Arrival Day (MJD), Angular Error, RA, Dec]

function create_sample_nus(orig_data::Array{Float64,2},N_nus=0)

    if N_nus == 0
        N_nus = length(orig_data[:,1])
    end

    mjd = sample(orig_data[:,1],replace=true,N_nus);
    ang_err = sample(orig_data[:,2],replace=true,N_nus);
    ra = rand(Uniform(0,2pi),N_nus);
    dec = sample(orig_data[:,4],replace=true,N_nus);

    return hcat(mjd,ang_err,ra,dec)
end

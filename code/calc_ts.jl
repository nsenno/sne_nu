using StatsBase
using Distributions
using Optim

code_dir = string(homedir(),"/DataScience/sne_nu/code/");
data_dir = string(homedir(),"/DataScience/sne_nu/data/");
notebook_dir = string(homedir(),"/DataScience/sne_nu/notebook/");

include(string(code_dir,"sig_and_bkg_pdfs.jl"));

# Function that calculates the TS for an array of s_j's and their associated
# sn objects
#
# Inputs :
#   nns  --  an array of s_j's, they do not need to maximize TS
#   t_sn --  an array of sn objects that correspond to the s_j array
#
# Output :
#   the value of TS corresponding to the s_j's and coefficients contained in
#   the array of sn objects

function calc_TS(nns::Array{Float64,1},t_sn::Array{sn,1})
		@assert length(nns) == length(t_sn)
    inner_sum = [sum(log(nns[i].*t_sn[i].coefs+ 1.0)) for i in 1:length(t_sn)]

    return sum( nns .- inner_sum);
end

# Function that calculates the derivative of TS. This is used to find the
# stationary points, which ultimately locate the maximum of the original TS.
# It turns out that the derivatives are separable, and can be calculated for
# each SN individually.
#
# Inputs :
#   t_n  --  the independent variable s_j for a single SN
#   t_cs --  an array of coefficients corresponding to a single SN
#
# Outputs :
#   the derivative of TS for a given SN
#

function calc_TS_deriv(t_n::Float64, t_cs::Array{Float64,1})
    return sum(t_cs./(t_n.*t_cs+1.0))-1.0
end

nu_data = readdlm(string(data_dir,"cleaned_nu_data.csv"),',',header=true)[1];
sne_data = readdlm(string(data_dir,"cleaned_sne_data.csv"),',',header=true)[1];

nus = create_nu_array(nu_data);
sne = create_sn_array(sne_data);

ns = zeros(Float64,length(sne));

nbs = readdlm(string(data_dir,"nb_data.dat"))[:];

[add_nb!(sn,nbs[i]) for (i,sn) in enumerate(sne)];

[add_nu!(sn,nus[find_associated_nus(sn,nus)]) for (i,sn) in enumerate(sne)];

map(calc_coefs!,sne);

# ensure that the sn that we are considering have associated nus, or else
# the calculation won't work

non_empty_arrays = [length(sn.coefs) != 0.0 for sn in sne]

# find the s_j's that form a stationary solution (i.e., the derivative of TS = 0)

ns[non_empty_arrays] = [optimize(x-> abs(calc_TS_deriv(x,sn.coefs)),max(maximum(-1./sn.coefs),-sn.nb),100.0).minimizer for sn in sne[non_empty_arrays]];

# calculate the value of TS corresponding to this array of s_j's

result_TS = calc_TS(ns[non_empty_arrays],sne[non_empty_arrays]);

println(result_TS)

using StatsBase
using Distributions
using Optim

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

# Function that returns the number of signal neutrinos for an array of sn
# objects, and the corresponding value of the TS
#
# Input :
#   sne  --  an array of sn objects that has been loaded with associated
#            nus and coefficients
#
# Outputs :
#   ns        --  array of s_j's corresponding to the input sne
#   result_TS --  resulting value of TS from optimization 

function get_TS(sne::Array{sn,1})
  # ensure that the sn that we are considering have associated nus, or else
  # the calculation won't work

  non_empty_arrays = [length(sn.coefs) != 0.0 for sn in sne]

  # find the s_j's that form a stationary solution (i.e., the derivative of TS = 0)

  ns[non_empty_arrays] = [optimize(x-> abs(calc_TS_deriv(x,sn.coefs)),max(maximum(-1./sn.coefs),-sn.nb),100.0).minimizer for sn in sne[non_empty_arrays]];

  # calculate the value of TS corresponding to this array of s_j's

  result_TS = calc_TS(ns[non_empty_arrays],sne[non_empty_arrays]);

  return ns, result_TS
end

# Code that finds all neutrinos not associated with SNe, and creates a
# PDF using Kernel Density Esitmation. This will serve as part of the
# background PDF for neutrino arrival direction

using Interpolations
using StatsBase
using KernelDensity

include("sn.jl")

# Function that determines if a neutrino is associated with any sn
#
# Input :
#   nu_idx  --  numerical index corresponding to a single neutrino
#   sn_list --  a collection of neutrino indeces, each corresponding to
#               the neutrinos associated with each sn
# Output :
#   Bool determining if neutrino is associated with any sn

function is_nu_in_dict(nu_idx, sn_list)
    for a_sn in sn_list
        if nu_idx in a_sn
            return true
        end
    end

    return false
end

nu_data = readdlm("data/cleaned_nu_data.csv",',',header=true)
sne_data = readdlm("data/cleaned_sne_data.csv",',',header=true)

len_nus = length(nu_data[1][:,1])
len_sne = length(sne_data[1][:,1])

# create a range of indceies that reference all neutrinos
# this will make finding the associated neutrinos easier and faster
nu_idxs = range(1,len_nus)

# import all of the neutrino data into neutrino structures
nus = [nu(nu_data[1][i,:]...) for i in 1:len_nus];

# import all of the SNe data into sn structures
sne = [sn(sne_data[1][i,2:5]...) for i in 1:len_sne];

# create a dictionary of neutrinos associated with each sn. The keys
# are the sn names and the values are a collection of indecies determined by
# the find_associated_nus() function

my_dict = Dict(sne_data[1][i,1] => nu_idxs[find_associated_nus(a_sn,nus)] for (i,a_sn) in enumerate(sne))

# Determine an array of Bools indicating if a neutrino is a background event
# i.e., it is not associated with any sn

idxs_of_non_associated_nus = [~is_contained(idx,values(my_dict)) for idx in nu_idxs];

# Generate an array of dec's for background neutrinos

bkg_dec = [particle.dec for particle in nus[idxs_of_non_associated_nus]];

writedlm("data/bkg_nu_dec_values.dat",bkg_dec)

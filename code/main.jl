code_dir = string(homedir(),"/DataScience/sne_nu/code/");
data_dir = string(homedir(),"/DataScience/sne_nu/data/");
notebook_dir = string(homedir(),"/DataScience/sne_nu/notebook/");

include(string(code_dir,"nu.jl"));
include(string(code_dir,"sn.jl"));
include(string(code_dir,"sig_and_bkg_pdfs.jl"));
include(string(code_dir,"calc_ts.jl"));
include(string(code_dir,"analysis_TS.jl"))

nu_data = readdlm(string(data_dir,"cleaned_nu_data.csv"),',',header=true)[1];
sne_data = readdlm(string(data_dir,"cleaned_sne_data.csv"),',',header=true)[1];

sne = create_sn_array(sne_data);

nbs = readdlm(string(data_dir,"nb_data.dat"))[:];

[add_nb!(sn,nbs[i]) for (i,sn) in enumerate(sne)];

load_sne!(nu_data,sne)

TS_obs = get_TS(sne)[2]

#println(get_TS(sne)[2])

#println(get_sample_TS(sne)[2])

sample_TSs = get_many_sample_TS(5);
println(sample_TSs)

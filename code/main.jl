code_dir = string(homedir(),"/DataScience/sne_nu/code/");
data_dir = string(homedir(),"/DataScience/sne_nu/data/");
notebook_dir = string(homedir(),"/DataScience/sne_nu/notebook/");

include(string(code_dir,"nu.jl"));
include(string(code_dir,"sn.jl"));
include(string(code_dir,"sig_and_bkg_pdfs.jl"));
include(string(code_dir,"calc_ts.jl"));

nu_data = readdlm(string(data_dir,"cleaned_nu_data.csv"),',',header=true)[1];
sne_data = readdlm(string(data_dir,"cleaned_sne_data.csv"),',',header=true)[1];

nus = create_nu_array(nu_data);
sne = create_sn_array(sne_data);

nbs = readdlm(string(data_dir,"nb_data.dat"))[:];

[add_nb!(sn,nbs[i]) for (i,sn) in enumerate(sne)];

[add_nu!(sn,nus[find_associated_nus(sn,nus)]) for (i,sn) in enumerate(sne)];

map(calc_coefs!,sne);

println(get_TS(sne)[2])

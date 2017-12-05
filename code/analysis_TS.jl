function get_sample_TS(temp_sne::Array{sn,1})

  sample_nu_data = create_sample_nus(nu_data)

  load_sne!(create_sample_nus(nu_data),temp_sne)

  sample_s_j, sample_TS = get_TS(temp_sne)

  return sample_s_j, sample_TS

end

function get_many_sample_TS(n_samples::Int, return_s_j = false)
  @assert n_samples > 0
  if !return_s_j
    return [get_sample_TS(sne)[2] for n in 1:n_samples]
  else
    return [get_sample_TS(sne) for n in 1:n_samples]
  end
end

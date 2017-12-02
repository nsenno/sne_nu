# Create an immutable type nu which contains the arrival time and direction (Right Assension RA and Declination Dec)
# These objects will be useful when creating random samples 

immutable nu
    mjd::Float64
    ang_err::Float64
    ra::Float64
    dec::Float64
end

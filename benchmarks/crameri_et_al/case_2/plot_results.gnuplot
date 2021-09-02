# Plot results from Crameri et. al. case 1 benchmark
# Includes curves for Underworld, SULEC, MILAMIN_VEP, as well
# as an approximate analytical solution

years_in_myrs = 1e-6

plot "UNDERWORLD_fs.dat" using 1:2, \
     "MILAMIN_fs.dat" using ($1/1000):2, \
     "SULEC_fs.dat" using 1:2, \
     "output/statistics.dat" using (column("Time (years)")*years_in_myrs):(column("Maximum topography (m)"))

pause -1


# Plot results from Crameri et. al. case 1 benchmark
# Includes curves for Underworld, SULEC, MILAMIN_VEP, as well
# as an approximate analytical solution

years_in_kyrs = 0.001
m_in_km = 0.001

plot "UNDERWORLD_fs.dat" using 1:2, \
     "MILAMIN_VEP.dat" using 1:2, \
     "SULEC_fs.dat" using 1:2, \
     "output/statistics.dat" using (column("Time (years)")*years_in_kyrs):(column("Maximum topography (m)")*m_in_km), \
     "output/statistics.dat" using (column("Time (years)")*years_in_kyrs):(7. * exp(column("Time (years)") * -1/14825)) with lines

pause -1


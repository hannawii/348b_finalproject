import sys
import math

color = sys.argv[6]
hemoglobin_frac = float(sys.argv[1])
melanin_frac = float(sys.argv[2])
melanin_blend = float(sys.argv[3])
oiliness = float(sys.argv[4])
gamma = float(sys.argv[5])
wavelength = 650

if (color == "red" or color == "r") : 
	wavelength = 650
	sigma_oxy = 0.5
	sigma_deoxy = 3.5
elif (color == "green" or color == "g") : 
	wavelength = 510
	sigma_oxy = 10
	sigma_deoxy = 20
elif (color == "blue" or color == "b") : 
	wavelength = 475
	sigma_oxy = 10.2
	sigma_deoxy = 9.5
else : print "Improper wavelength, default to red."

sigma_em = (6.6 * math.pow(10, 10)) * math.pow(wavelength, -3.33)
sigma_pm = (2.9 * math.pow(10, 14)) * math.pow(wavelength, -4.75)
sigma_base = 0.0244 + (8.53 * math.exp(-(wavelength - 154) / 66.2))
sigma_epi = melanin_frac * (melanin_blend * sigma_em + (1 - melanin_blend) * sigma_pm) + (1 - melanin_frac) * sigma_base

sigma_derm = hemoglobin_frac * (gamma * sigma_oxy + (1 - gamma) * sigma_deoxy) + (1 - hemoglobin_frac) * sigma_base

sigma_prime_s_epi = 14.74 * math.pow(wavelength, -0.22) + (2.2 * math.pow(10, 11) * math.pow(wavelength, -4))
sigma_prime_s_derm = sigma_prime_s_epi * 0.5

print "sigma_a_epi", sigma_epi, "\nsigma_a_derm", sigma_derm, "\nsigma_prime_s_epi", sigma_prime_s_epi, "\nsigma_prime_s_derm", sigma_prime_s_derm




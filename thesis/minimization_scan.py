# Program to make a scan over many minimizations, and investigate which give the lowest best-fit point
import numpy as np 


# Give starting points as lists
sp_sq_list = []
sp_chi2_list = []
sp_sl_list = []
sp_chi1_list = []

# Loop over all starting point combinations, omitting unphysical ones if there are any
for sp_sq in sp_sq_list:

	for sp_chi2 in sp_chi2_list:
		if sp_chi2 >= sp_sq:
			continue

		for sp_sl in sp_sl_list:
			if sp_sl >= sp_chi2:
				continue

			for sp_chi1 in sp_chi1_list:
				if sp_chi1 >= sp_sl:
					continue

				# RUN MINIMIZATION
				
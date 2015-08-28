#! /usr/bin/env python
import os

def rewrite_species_result(input_species_file, output_species_file):
	try:
		with open(input_species_file) as f:
			content = f.read().splitlines()
		f.close()

		largestSpecies = 0
		taxaList = []
		assignments = {}

		for i in range(1,151):
			line = content[i]
			line = " ".join(line.split())
			species_idx = int(line.split(' ')[1])
			taxon_name = line.split(' ')[2]
			if (species_idx > largestSpecies):
				assignments[species_idx] = []
				largestSpecies = species_idx
			assignments[species_idx].append(taxon_name)

		if not os.path.exists(os.path.dirname(output_species_file)):
    			os.makedirs(os.path.dirname(output_species_file))
		speciesOut = open(output_species_file, 'w')

		speciesOut.write("Species 1:\n")
		for j in range(0, len(assignments[1])):
				speciesOut.write(assignments[1][j] + "\n")
		for i in range(2, largestSpecies + 1):
			speciesOut.write("\nSpecies " + str(i) + ":\n")
			for j in range(0, len(assignments[i])):
				speciesOut.write(assignments[i][j] + "\n")

		speciesOut.close()
	except IOError:
		print "File not found: " + input_species_file

set_names = ["Ne10000", "Ne100000", "Ne500000", "Ne1000000"]

#rewrite_species_result("gmyc_results/set_Ne10000/gmyc_results_set_Ne10000.1.txt", "similar_to_GMYC_gmyc/set_Ne10000/gmyc_results_set_Ne10000.1.txt")

for set_name in set_names:
	for i in range(1,101):
		input_species_file = "gmyc_results/set_" + set_name + "/gmyc_results_set_" + set_name + "." + str(i) + ".txt"
		output_species_file = "similar_to_GMYC_gmyc_minbr_0/set_" + set_name + "/gmyc_results_set_" + set_name + "." + str(i) + ".txt"
		rewrite_species_result(input_species_file, output_species_file)

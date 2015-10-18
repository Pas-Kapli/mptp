#! /usr/bin/env python
import os

def rewrite_species_result(input_species_file, output_species_file):
	try:
		with open(input_species_file) as f:
			content = f.read().splitlines()
		f.close()

		taxaListString = content[0].split(':')[1]
		taxaList = taxaListString.split(',')
		speciesList = content[1].split(',')

		if not os.path.exists(os.path.dirname(output_species_file)):
    			os.makedirs(os.path.dirname(output_species_file))
		speciesOut = open(output_species_file, 'w')

		oldSpeciesIdx = speciesList[0]
		speciesOut.write("Species " + oldSpeciesIdx + ":\n")
		speciesOut.write(taxaList[0] + "\n")
		for i in range(1,len(speciesList)):
			if (speciesList[i] == oldSpeciesIdx):
				speciesOut.write(taxaList[i] + "\n")
			else:
				oldSpeciesIdx = speciesList[i]
				speciesOut.write("\nSpecies " + oldSpeciesIdx + ":\n")
				speciesOut.write(taxaList[i] + "\n")

		speciesOut.close()
	except IOError:
		print "File not found: " + input_species_file

set_names = ["Ne10000", "Ne100000", "Ne500000", "Ne1000000"]

for set_name in set_names:
	for i in range(1,101):
		if set_name == "Ne10000":
			input_species_file = "similar_to_GMYC/15-08-2015.16-40/set_BIRTH0.27_" + set_name + "/PTP_result_BIRTH0.27_" + set_name + "_" + str(i) + ".PTPPartitions.txt"
		elif set_name == "Ne500000":
			input_species_file = "similar_to_GMYC/15-08-2015.16-40/set_BIRTH0.27_" + set_name + "/PTP_BIRTH0.27_" + set_name + "_" + str(i) + ".PTPPartitions.txt"
		elif set_name == "Ne100000":
			input_species_file = "similar_to_GMYC/15-08-2015.16-40/set_BIRTH0.27_" + set_name + "/PTP_BIRTH0.27_" + set_name + "_" + str(i) + ".PTPPartitions.txt"
		else:
			input_species_file = "similar_to_GMYC/15-08-2015.16-40/set_BIRTH0.27_" + set_name + "/PTP_result." + str(i) + ".PTPPartitions.txt"
		output_species_file = "similar_to_GMYC_PTP_minbr_default/set_" + set_name + "/PTP_results_set_" + set_name + "." + str(i) + ".txt"
		rewrite_species_result(input_species_file, output_species_file)

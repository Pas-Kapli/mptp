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

set_names = ["1", "5", "10", "20", "40", "80", "160"]

for set_name in set_names:
	for i in range(1,100):
		if (set_name == "1"):
			input_species_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/PTP_results_set" + set_name + "/PTP_result_" + set_name + "." + str(i) + ".PTPPartitions.txt"
		elif (set_name == "80"):
			input_species_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/PTP_results_set" + set_name + "_unique_taxa/PTP_result_rerooted.unique_taxa_" + set_name + "." + str(i) + ".PTPPartitions.txt"
		elif (set_name == "40"):
			input_species_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/PTP_results_set" + set_name + "/PTP_result_rerooted.unique_taxa_" + set_name + "." + str(i) + ".PTPPartitions.txt"
		elif (set_name == "5"):
			input_species_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/PTP_results_set" + set_name + "/PTP_result_rerooted.unique_taxa_" + set_name + "." + str(i) + ".PTPPartitions.txt"
		else:
			input_species_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/PTP_results_set" + set_name + "/PTP_result_rerooted.unique_taxa." + set_name + "." + str(i) + ".PTPPartitions.txt"
		output_species_file = "unique_taxa_big_PTP/set_" + set_name + "/PTP_results_set_" + set_name + "." + str(i) + ".txt"
		rewrite_species_result(input_species_file, output_species_file)

#! /usr/bin/env python
import os

def create_subsets(alignmentFile, num_of_species, sum_of_species, num_basepairs, output_taxa_file, output_alignment_file, num_alignments):
	try:
		with open(alignmentFile) as f:
			content = f.read().splitlines()
		f.close()
		speciesList = []

		for i in range(0,31):
			emptyList = []
			speciesList.append(emptyList)

		alignments = {}
		for i in range(1, len(content)): # ignore first line
			contentSplitted = content[i].split();
			taxonName = contentSplitted[0]
			alignments[taxonName] = contentSplitted[1][0:num_basepairs]
			species = taxonName.split('.')[0]
			speciesList[int(species)].append(taxonName)

		speciesListSorted = sorted(speciesList, key = len)

		currentIdx = 0

		selectedTaxa = []

		found = 0
		for i in range(30,-1,-1):
			if currentIdx < len(num_of_species):
				if len(speciesListSorted[i]) >= sum_of_species[currentIdx]:
					found = found + 1
					for j in range(1, sum_of_species[currentIdx]):
						selectedTaxa.append(speciesListSorted[i][j])
				else:
					print "We had an error :("
				if found == num_of_species[currentIdx]:
					currentIdx = currentIdx + 1
					found = 0

		# write the solutions into the files
		if not os.path.exists(os.path.dirname(output_taxa_file)):
    			os.makedirs(os.path.dirname(output_taxa_file))
		taxaOut = open(output_taxa_file, 'w')
		for taxon in selectedTaxa:
			taxaOut.write(taxon + "\n")
		taxaOut.close()

		if not os.path.exists(os.path.dirname(output_alignment_file)):
    			os.makedirs(os.path.dirname(output_alignment_file))
		alignmentOut = open(output_alignment_file, 'w')
		alignmentOut.write(str(num_alignments) + " " + str(num_basepairs) + "\n")
		for taxon in selectedTaxa:
			alignmentOut.write(taxon + " "+ alignments[taxon] + "\n")
		alignmentOut.close()

		return (currentIdx >= len(num_of_species))
    	except IOError:
		print "File not found: " + alignmentFile


set_names = ["set_1", "set_5", "set_10", "set_20", "set_40", "set_80", "set_160"]
num_of_species = [3, 6, 9, 12]
size_of_species = [35, 25, 10, 2]
uniform_num = [30]
uniform_size = [12]
base_pairs = [100, 250, 500, 1000]
uniform_num_alignments = 360
nonuniform_num_alignments = 369

for set_name in set_names:
	for i in range(1,100):
		for bp in base_pairs:
			output_nonuniform_taxa_file = "nonuniform/taxa/"+str(bp)+"/taxa.simulated_" + set_name + "_" + str(i)
			output_nonuniform_alignment_file = "nonuniform/alignments/"+str(bp)+"/simulated_tree_" + set_name + "_" + str(i)
			output_uniform_taxa_file = "uniform/taxa/"+str(bp)+"/taxa.simulated_" + set_name + "_" + str(i)
			output_uniform_alignment_file = "uniform/alignments/"+str(bp)+"/simulated_tree_" + set_name + "_" + str(i)

			alignmentFile = "reduced_alignments/" + set_name + "/simulated_" + set_name + "_" + str(i) + ".phy.reduced"
			if create_subsets(alignmentFile, num_of_species, size_of_species, bp, output_nonuniform_taxa_file, output_nonuniform_alignment_file, nonuniform_num_alignments) == False:
				print "Found a file that does not fit our requirement :-("
			if create_subsets(alignmentFile, uniform_num, uniform_size, bp, output_uniform_taxa_file, output_uniform_alignment_file, uniform_num_alignments) == False:
				print "Found a file that does not fit our requirement :-("

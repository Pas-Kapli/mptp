#! /usr/bin/env python
import os
import commands

def extract_tree_score(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Tree penalty score:"):
			return int(line.split(': ')[1])
			break

def extract_nmi_score(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("NMI score:"):
			return float(line.split(': ')[1])
			break

def extract_num_species(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Number of species in input file:"):
			return int(line.split(': ')[1])
			if (int(line.split(': ')[1]) == 1):
				print "Baaaaad data"
			break

def extract_num_real_species(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Number of real species:"):
			return int(line.split(': ')[1])
			break

def extract_score_real_single(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Score real single:"):
			return float(line.split(': ')[1])
			break

def extract_score_real_multi(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Score real multi:"):
			return float(line.split(': ')[1])
			break

def extract_score_input_single(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Score input single:"):
			return float(line.split(': ')[1])
			break

def extract_score_input_multi(input_text):
	lines = input_text.split('\n')
	for line in lines:
		if line.startswith("Score input multi:"):
			return float(line.split(': ')[1])
			break

def grab_scorings(input_tree_file, output_delimit_single_minbr_0, output_delimit_multi_minbr_0, output_delimit_single_minbr_default, output_delimit_multi_minbr_default, output_PTP_minbr_0):
	try:
		open(input_tree_file)
		programNames = ['delimit_single_minbr_0', 'delimit_multi_minbr_0', 'delimit_single_minbr_default', 'delimit_multi_minbr_default', 'PTP_minbr_0']
		tree_scores = {}
		nmi_scores = {}
		num_species = {}
		single_scores = {}
		multi_scores = {}
		num_real_species = 0
		score_real_single = 0
		score_real_multi = 0

		tree_scores['delimit_single_minbr_0'] = extract_tree_score(output_delimit_single_minbr_0)
		tree_scores['delimit_multi_minbr_0'] = extract_tree_score(output_delimit_multi_minbr_0)
		tree_scores['delimit_single_minbr_default'] = extract_tree_score(output_delimit_single_minbr_default)
		tree_scores['delimit_multi_minbr_default'] = extract_tree_score(output_delimit_multi_minbr_default)
		tree_scores['PTP_minbr_0'] = extract_tree_score(output_PTP_minbr_0)

		nmi_scores['delimit_single_minbr_0'] = extract_nmi_score(output_delimit_single_minbr_0)
		nmi_scores['delimit_multi_minbr_0'] = extract_nmi_score(output_delimit_multi_minbr_0)
		nmi_scores['delimit_single_minbr_default'] = extract_nmi_score(output_delimit_single_minbr_default)
		nmi_scores['delimit_multi_minbr_default'] = extract_nmi_score(output_delimit_multi_minbr_default)
		nmi_scores['PTP_minbr_0'] = extract_nmi_score(output_PTP_minbr_0)

		num_species['delimit_single_minbr_0'] = extract_num_species(output_delimit_single_minbr_0)
		num_species['delimit_multi_minbr_0'] = extract_num_species(output_delimit_multi_minbr_0)
		num_species['delimit_single_minbr_default'] = extract_num_species(output_delimit_single_minbr_default)
		num_species['delimit_multi_minbr_default'] = extract_num_species(output_delimit_multi_minbr_default)
		num_species['PTP_minbr_0'] = extract_num_species(output_PTP_minbr_0)

		single_scores['delimit_single_minbr_0'] = extract_score_input_single(output_delimit_single_minbr_0)
		single_scores['delimit_multi_minbr_0'] = extract_score_input_single(output_delimit_multi_minbr_0)
		single_scores['delimit_single_minbr_default'] = extract_score_input_single(output_delimit_single_minbr_default)
		single_scores['delimit_multi_minbr_default'] = extract_score_input_single(output_delimit_multi_minbr_default)
		single_scores['PTP_minbr_0'] = extract_score_input_single(output_PTP_minbr_0)

		multi_scores['delimit_single_minbr_0'] = extract_score_input_multi(output_delimit_single_minbr_0)
		multi_scores['delimit_multi_minbr_0'] = extract_score_input_multi(output_delimit_multi_minbr_0)
		multi_scores['delimit_single_minbr_default'] = extract_score_input_multi(output_delimit_single_minbr_default)
		multi_scores['delimit_multi_minbr_default'] = extract_score_input_multi(output_delimit_multi_minbr_default)
		multi_scores['PTP_minbr_0'] = extract_score_input_multi(output_PTP_minbr_0)

		score_real_single = extract_score_real_single(output_delimit_single_minbr_0)
		score_real_multi = extract_score_real_multi(output_delimit_single_minbr_0)
		num_real_species = extract_num_real_species(output_delimit_single_minbr_0)

		return (tree_scores, nmi_scores, num_species, single_scores, multi_scores, score_real_single, score_real_multi, num_real_species)
	except IOError:
		print "File not found: " + input_tree_file

def create_scoring_results(input_tree_file, input_delimit_single_minbr_0_file, input_delimit_multi_minbr_0_file, input_delimit_single_minbr_default_file, input_delimit_multi_minbr_default_file, input_PTP_minbr_0_file, output_delimit_single_minbr_0_file, output_delimit_multi_minbr_0_file, output_delimit_single_minbr_default_file, output_delimit_multi_minbr_default_file, output_PTP_minbr_0_file):
	try:
		open(input_tree_file)
		if not os.path.exists(os.path.dirname(output_delimit_single_minbr_0_file)):
	    			os.makedirs(os.path.dirname(output_delimit_single_minbr_0_file))
		if not os.path.exists(os.path.dirname(output_delimit_multi_minbr_0_file)):
	    			os.makedirs(os.path.dirname(output_delimit_multi_minbr_0_file))
		if not os.path.exists(os.path.dirname(output_delimit_single_minbr_default_file)):
	    			os.makedirs(os.path.dirname(output_delimit_single_minbr_default_file))
		if not os.path.exists(os.path.dirname(output_delimit_multi_minbr_default_file)):
	    			os.makedirs(os.path.dirname(output_delimit_multi_minbr_default_file))
		if not os.path.exists(os.path.dirname(output_PTP_minbr_0_file)):
	    			os.makedirs(os.path.dirname(output_PTP_minbr_0_file))

		call_delimit_single_minbr_0 = "./delimit --score " + input_delimit_single_minbr_0_file + " --min_br 0 --tree_file " + input_tree_file + " --output_file foo"
		call_delimit_multi_minbr_0 = "./delimit --score " + input_delimit_multi_minbr_0_file + " --min_br 0 --tree_file " + input_tree_file + " --output_file foo"
		call_delimit_single_minbr_default = "./delimit --score " + input_delimit_single_minbr_default_file + " --tree_file " + input_tree_file + " --output_file foo"
		call_delimit_multi_minbr_default = "./delimit --score " + input_delimit_multi_minbr_default_file + " --tree_file " + input_tree_file + " --output_file foo"
		call_PTP_minbr_0 = "./delimit --score " + input_PTP_minbr_0_file + " --min_br 0 --tree_file " + input_tree_file + " --output_file foo"

		(stat_delimit_single_minbr_0, output_delimit_single_minbr_0) = commands.getstatusoutput(call_delimit_single_minbr_0)
		(stat_delimit_multi_minbr_0, output_delimit_multi_minbr_0) = commands.getstatusoutput(call_delimit_multi_minbr_0)
		(stat_delimit_single_minbr_default, output_delimit_single_minbr_default) = commands.getstatusoutput(call_delimit_single_minbr_default)
		(stat_delimit_multi_minbr_default, output_delimit_multi_minbr_default) = commands.getstatusoutput(call_delimit_multi_minbr_default)
		(stat_PTP_minbr_0, output_PTP_minbr_0) = commands.getstatusoutput(call_PTP_minbr_0)

		delimit_single_minbr_0_out = open(output_delimit_single_minbr_0_file, 'w')
		delimit_multi_minbr_0_out = open(output_delimit_multi_minbr_0_file, 'w')
		delimit_single_minbr_default_out = open(output_delimit_single_minbr_default_file, 'w')
		delimit_multi_minbr_default_out = open(output_delimit_multi_minbr_default_file, 'w')
		PTP_minbr_0_out = open(output_PTP_minbr_0_file, 'w')

		delimit_single_minbr_0_out.write(output_delimit_single_minbr_0)
		delimit_multi_minbr_0_out.write(output_delimit_multi_minbr_0)
		delimit_single_minbr_default_out.write(output_delimit_single_minbr_default)
		delimit_multi_minbr_default_out.write(output_delimit_multi_minbr_default)
		PTP_minbr_0_out.write(output_PTP_minbr_0)

		delimit_single_minbr_0_out.close()
		delimit_multi_minbr_0_out.close()
		delimit_single_minbr_default_out.close()
		delimit_multi_minbr_default_out.close()
		PTP_minbr_0_out.close()

		return grab_scorings(input_tree_file, output_delimit_single_minbr_0, output_delimit_multi_minbr_0, output_delimit_single_minbr_default, output_delimit_multi_minbr_default, output_PTP_minbr_0)
	except IOError:
		print "File not found: " + input_tree_file

set_names = ["1", "5", "10", "20", "40", "80", "160"]
names = ['delimit_single_minbr_0', 'delimit_multi_minbr_0', 'delimit_single_minbr_default', 'delimit_multi_minbr_default', 'PTP_minbr_0']

gnuplotOut_tree_scores = open('workfile_tree_scores', 'w')
gnuplotOut_nmi_scores = open('workfile_nmi_scores', 'w')
gnuplotOut_single_scores = open('workfile_single_scores', 'w')
gnuplotOut_multi_scores = open('workfile_multi_scores', 'w')
gnuplotOut_num_species = open('workfile_num_species', 'w')

for set_name in set_names:
	num_valid_indices = 0
	average_tree_scores = {}	
	average_nmi_scores = {}
	average_num_species = {}
	average_single_scores = {}
	average_multi_scores = {}
	average_real_num_species = 0
	average_real_score_single = 0
	average_real_score_multi = 0

	for name in names:
		average_tree_scores[name] = 0
		average_nmi_scores[name] = 0
		average_num_species[name] = 0
		average_single_scores[name] = 0
		average_multi_scores[name] = 0

	for i in range(1,100):
		if (set_name == "1"):
			input_tree_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/rooted.inferred_unique_taxa." + str(i)
		else:
			input_tree_file = "unique_taxa_trees_big_dataset/set_" + set_name + "/RAxML_inferred_trees_unique_taxa/rooted.inferred_unique_taxa_set_" + set_name + "." + str(i)

		try:
			open(input_tree_file)

			input_delimit_single_minbr_0_file = "unique_taxa_big_delimit_single_minbr_0/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
			input_delimit_multi_minbr_0_file = "unique_taxa_big_delimit_multi_minbr_0/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
			input_delimit_single_minbr_default_file = "unique_taxa_big_delimit_single_minbr_default/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
			input_delimit_multi_minbr_default_file = "unique_taxa_big_delimit_multi_minbr_default/set_" + set_name + "/delimit_results_set_" + set_name + "." + str(i) + ".txt"
			input_PTP_minbr_0_file = "unique_taxa_big_PTP_minbr_0/set_" + set_name + "/PTP_results_set_" + set_name + "." + str(i) + ".txt"
		
			score_path = "unique_taxa_big_scoring_results/"		
			output_delimit_single_minbr_0_file = score_path + "delimit_single_minbr_0/set_" + set_name + "/delimit_score_set_" + set_name + "." + str(i) + ".txt"
			output_delimit_multi_minbr_0_file = score_path + "delimit_multi_minbr_0/set_" + set_name + "/delimit_score_set_" + set_name + "." + str(i) + ".txt"
			output_delimit_single_minbr_default_file = score_path + "delimit_single_minbr_default/set_" + set_name + "/delimit_score_set_" + set_name + "." + str(i) + ".txt"
			output_delimit_multi_minbr_default_file = score_path + "delimit_multi_minbr_default/set_" + set_name + "/delimit_score_set_" + set_name + "." + str(i) + ".txt"
			output_PTP_minbr_0_file = score_path + "PTP_minbr_0/set_" + set_name + "/PTP_score_set_" + set_name + "." + str(i) + ".txt"

			(tree_scores, nmi_scores, num_species, single_scores, multi_scores, score_real_single, score_real_multi, num_real_species) = create_scoring_results(input_tree_file, input_delimit_single_minbr_0_file, input_delimit_multi_minbr_0_file, input_delimit_single_minbr_default_file, input_delimit_multi_minbr_default_file, input_PTP_minbr_0_file, output_delimit_single_minbr_0_file, output_delimit_multi_minbr_0_file, output_delimit_single_minbr_default_file, output_delimit_multi_minbr_default_file, output_PTP_minbr_0_file)

			try:
				for name in names:
					average_tree_scores[name] = average_tree_scores[name] + tree_scores[name]
					average_nmi_scores[name] = average_nmi_scores[name] + nmi_scores[name]
					average_num_species[name] = average_num_species[name] + num_species[name]
					average_single_scores[name] = average_single_scores[name] + single_scores[name]
					average_multi_scores[name] = average_multi_scores[name] + multi_scores[name]
				average_real_num_species = average_real_num_species + num_real_species
				average_real_score_single = average_real_score_single + score_real_single
				average_real_score_multi = average_real_score_multi + score_real_multi
			except:
				print "File is empty: " + input_tree_file
				num_valid_indices = num_valid_indices - 1

			num_valid_indices = num_valid_indices + 1
		except IOError:
			#1
			print "File not found: " + input_tree_file
	for name in names:
		average_tree_scores[name] = float(average_tree_scores[name]) / float(num_valid_indices)
		average_nmi_scores[name] = float(average_nmi_scores[name]) / float(num_valid_indices)
		average_num_species[name] = float(average_num_species[name]) / float(num_valid_indices)
		average_single_scores[name] = float(average_single_scores[name]) / float(num_valid_indices)
		average_multi_scores[name] = float(average_multi_scores[name]) / float(num_valid_indices)
		
		#print "Set " + set_name + ": Average tree score " + name
		#print average_tree_scores[name]
		#print "Set " + set_name + ": Average NMI score " + name
		#print average_nmi_scores[name]
		#print "Set " + set_name + ": Average num species " + name
		#print average_num_species[name]
		#print "Set " + set_name + ": Average input score single " + name
		#print average_single_scores[name]
		#print "Set " + set_name + ": Average input score multi " + name
		#print average_multi_scores[name]
	average_real_num_species = float(average_real_num_species) / float(num_valid_indices)
	average_real_score_single = float(average_real_score_single) / float(num_valid_indices)
	average_real_score_multi = float(average_real_score_multi) / float(num_valid_indices)
	#print "Set " + set_name + ": Average real num species "
	#print average_real_num_species
	#print "Set " + set_name + ": Average real score single "
	#print average_real_score_single
	#print "Set " + set_name + ": Average real score multi "
	#print average_real_score_multi

	gnuplotOut_tree_scores.write(set_name + ' ' + str(average_tree_scores['delimit_single_minbr_0']) + ' ' + str(average_tree_scores['delimit_multi_minbr_0']) + ' ' + str(average_tree_scores['delimit_single_minbr_default']) + ' ' + str(average_tree_scores['delimit_multi_minbr_default']) + ' ' + str(average_tree_scores['PTP_minbr_0']) + '\n')

	gnuplotOut_nmi_scores.write(set_name + ' ' + str(average_nmi_scores['delimit_single_minbr_0']) + ' ' + str(average_nmi_scores['delimit_multi_minbr_0']) + ' ' + str(average_nmi_scores['delimit_single_minbr_default']) + ' ' + str(average_nmi_scores['delimit_multi_minbr_default']) + ' ' + str(average_nmi_scores['PTP_minbr_0']) + '\n')

	gnuplotOut_single_scores.write(set_name + ' ' + str(average_single_scores['delimit_single_minbr_0']) + ' ' + str(average_single_scores['delimit_multi_minbr_0']) + ' ' + str(average_single_scores['delimit_single_minbr_default']) + ' ' + str(average_single_scores['delimit_multi_minbr_default']) + ' ' + str(average_single_scores['PTP_minbr_0']) + ' ' + str(average_real_score_single) + '\n')

	gnuplotOut_multi_scores.write(set_name + ' ' + str(average_multi_scores['delimit_single_minbr_0']) + ' ' + str(average_multi_scores['delimit_multi_minbr_0']) + ' ' + str(average_multi_scores['delimit_single_minbr_default']) + ' ' + str(average_multi_scores['delimit_multi_minbr_default']) + ' ' + str(average_multi_scores['PTP_minbr_0']) + ' ' + str(average_real_score_multi) + '\n')

	gnuplotOut_num_species.write(set_name + ' ' + str(average_num_species['delimit_single_minbr_0']) + ' ' + str(average_num_species['delimit_multi_minbr_0']) + ' ' + str(average_num_species['delimit_single_minbr_default']) + ' ' + str(average_num_species['delimit_multi_minbr_default']) + ' ' + str(average_num_species['PTP_minbr_0']) + ' ' + str(average_real_num_species) + '\n')

gnuplotOut_tree_scores.close()
gnuplotOut_nmi_scores.close()
gnuplotOut_single_scores.close()
gnuplotOut_multi_scores.close()
gnuplotOut_num_species.close()

commands.getstatusoutput('gnuplot plotscript')

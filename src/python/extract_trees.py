#! /usr/bin/env python
import os
import commands

set_names = ["Ne1e+05", "Ne1e+06", "Ne5e+05", "Ne10000"]

for set_name in set_names:
	try:
		tree_path = "SimulB&C." + set_name + "_nospec.phy"
		tree_file = open(tree_path)
		lines = tree_file.readlines()

		for i in range(1,101): # only the first 100 trees
			tree_destination = "SimulB&C_trees/set_" + set_name + "/SimulB&C_tree_set_" + set_name + "." + str(i) + ".txt"
			if not os.path.exists(os.path.dirname(tree_destination)):
		    		os.makedirs(os.path.dirname(tree_destination))
			tree_destination_file = open(tree_destination, 'w')
			print(lines[i-1])
		
		tree_file.close()
	except IOError:
		print "File not found: " + tree_path

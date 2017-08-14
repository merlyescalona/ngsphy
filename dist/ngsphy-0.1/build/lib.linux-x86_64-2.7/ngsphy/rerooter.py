#!/usr/bin/env python
# -*- coding: utf-8 -*-
import dendropy,logging,os,shutil,sys
import settings as sp

class Rerooter:
	"""
	Class used for the manipulation of gene trees, specifically, to handle
	the re-rooting and prunning of the given gene tree.
	----------------------------------------------------------------------------
	Attributes:
	- appLogger: logger to store status of the process flow
	- settings: Settings object withh all the program parameters
	- tree: final gene tree
	- outputFilename: filename where the new tree will be stored
	- outputFilePath: absolute path of the filename where the new tree will be stored
	"""
	# general
	appLoger=None
	settings=None
	# tree-related
	tree=None
	outputFilename="ngsphy.tree"
	outputFilePath=""

	def __init__(self, settings):
		self.appLogger=logging.getLogger('ngsphy')
		self.appLogger.debug('Rerooting')
		self.settings=settings
		self.outputFilePath=os.path.join(\
			self.settings.outputFolderPath,\
			self.outputFilename\
		)


	def run(self):
		"""
		Process flow of the re-rooting and pruninc process
		"""
		self.appLogger.debug('Running rerooting')
		try:
			self.tree=dendropy.Tree.get(path=self.settings.newickFilePath, schema="newick",preserve_underscores=True)
		except Exception as ex:
			return False, ex
		# print(self.tree.as_ascii_plot())
		newroot = self.tree.find_node_with_taxon_label(self.settings.referenceTipLabel)
		if newroot:
			self.tree.reroot_at_edge(newroot.edge, update_bipartitions=False)
			# print(self.tree.as_ascii_plot())
		else:
			return False, "{0}\n\t{1}\n\t{2}".format(\
			"Rerooting problem.",\
			"Something might be wrong with the reference label.",\
			"Please Verify. Exiting"\
			)

		leaves=[node.taxon.label for node in self.tree.leaf_node_iter() if not node.taxon.label == self.settings.referenceTipLabel]
		# print(leaves)
		try:
			mrca=self.tree.mrca(taxon_labels=leaves)
			self.tree.prune_taxa_with_labels([self.settings.referenceTipLabel])
			# It does not matter the naming of the  tips since this process will only
			# generate haploid individuals
			# for item in self.tree.leaf_node_iter():
			# 	if self.tree.seed_node==item.parent_node:
			# 		item.taxon.label="0_0_0"
			self.writeTreeIntoFile()
			self.settings.newickFilePath=self.outputFilePath
		except Exception as ex:
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="Rerooter - run: {0} | {1} - File: {2} - Line:{3}".format(ex,exc_type, fname, exc_tb.tb_lineno)
			return False, message

		return True,""

	def generateFolderStructure(self):
		"""
		Generation of basic folder structure for this process
		"""
		folder=os.path.join(self.settings.path,self.settings.projectName,"1")
		try:
			os.makedirs(folder)
			self.appLogger.info("Generating project folder ({0})".format(folder))
		except:
			self.appLogger.debug("Project folder exists ({0})".format(folder))
		try:
			shutil.copyfile(self.settings.referenceSequenceFilePath, os.path.joing(folder, "reference.fasta"))
		except:
			self.appLogger.debug("File already exists in this location ({0})".format(os.path.joing(folder, "reference.fasta")))

	def writeTreeIntoFile(self):
		"""
		Writes into a file the resulting tree
		"""
		self.tree.write(\
			path=self.outputFilePath,\
			schema="newick",\
			suppress_rooting=True\
			)

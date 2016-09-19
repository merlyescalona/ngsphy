#!/usr/bin/home/python

################################################################################
#
#    Copyright (C) 2016  Merly Escalona  <merlyescalona@uvigo.es>
#    SimPhyMating is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SimPhyMating is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=1
FIX_VERSION=0
PROGRAMNAME="simphy.ngs.wrapper.py"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
################################################################################
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
CONFIG_LEVEL=45
################################################################################
# python mating.py -p <prefixES_sequence_filename> -SF <simphy_path>
# python mating.py -p data -sf /media/merly/ME-Conus/conusSim/csTest1/
################################################################################
import argparse,datetime,logging,os,sys
import numpy as np
import random as rnd
from MEOutputFormatter import MEOutputFormatter as mof
from MELoggingFormatter import MELoggingFormatter as mlf
from Mating import Mating as snwm
from select import select


def handlingParameters():
    parser = argparse.ArgumentParser(\
        prog="{0} (v.{1}.{2}.{3})".format(PROGRAMNAME,VERSION,MIN_VERSION,FIX_VERSION),\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
        description=\
            '''
SimPhy NGS wrapper
==================

A program to simulate mating and next-generation sequencing reads from a SimPhy project.

For more detailed information, please check the README.md file.

            ''',\
        epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
        add_help=False
        )
    requiredGroup= parser.add_argument_group('Required arguments')
    requiredGroup.add_argument('-sf','--simphy-folder',metavar='<simphy_path>', type=str,\
        help='Path of the SimPhy project folder.', required=True)
    requiredGroup.add_argument('-p','--prefix', metavar='<dataset_prefix>', type=str,\
        help='Prefix of the filename that the sequences have', required=True)

    optionalGroup= parser.add_argument_group('Optional arguments')
    optionalGroup.add_argument('-nf','--number-files',metavar='<num_files>',type=int,\
        choices=range(1, 3),default=3,\
        help="Number of files generated per individual. 1: Only one file for both sequences. 2: A file per strand sequence. 3: A file with both sequences and a file per strand sequence. Default: 3")
    optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
        choices=LOG_LEVEL_CHOICES, default="INFO",\
        help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(LOG_LEVEL_CHOICES,LOG_LEVEL_CHOICES[1]))
    informationGroup= parser.add_argument_group('Information arguments')
    informationGroup.add_argument('-v', '--version',\
        action='version',\
        version='Mating version {0}.{1}.{2}'.format(VERSION,MIN_VERSION,FIX_VERSION),\
        help="Show program's version number and exit")
    informationGroup.add_argument('-h', '--help',\
        action='store_true',\
        help="Show this help message and exit")
    try:
        tmpArgs = parser.parse_args()
    except:
        sys.stdout.write("\n----------------------------------------------------------------------\n")
        parser.print_help()
        sys.exit()
    return tmpArgs

"""
###############################   MAIN   #######################################
"""
if __name__=="__main__":
    try:
        cmdArgs = handlingParameters()
        print(cmdArgs)
        prog = snwm.Mating(cmdArgs)
        prog.run()
    except KeyboardInterrupt:
        sys.stdout.write("{0}{1}\nInterrupted!{2}\nPlease run again for the expected outcome.\n".format(mof.BOLD,mof.DARKCYAN, mof.END))
        sys.exit()

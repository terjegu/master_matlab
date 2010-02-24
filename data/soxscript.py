#!/usr/bin/python
#Script for creating trainset.scp
"""create_trainset.py -- Make a list of wave files.

Switches:
    -h, --help         Displays this help message.
    -i <inputDir>,     --inputDir=<inputDir>
                       Input directory.
    -o <outputDir>,    --outputDir=<outputDir>
                       Output directory.

Examples:
    py soxscript.py -i source -o source_down """

import sys
import os
import getopt

def main():
	try:
		opts, xargs = getopt.getopt(sys.argv[1:], 'ho:i:o:', \
							['help', 'inputDir=', 'outputDir='])
	except getopt.error:
		print "Invalid switch!"
		sys.exit(1)
	# Print usage if no options are passed
	if opts == []:
		print __doc__
		sys.exit(1)

	for opt, arg in opts:
		if opt in ['-h', '--help']:
			print __doc__
			sys.exit(0)
		elif opt in ['-i','--inputDir']:
			inputDir = arg
		elif opt in ['-o','--outputDir']:
			outputDir = arg
	
	if os.path.isdir(inputDir) is False:
		print 'input directory does not exist'
		sys.exit(1)
		
	Makedir(outputDir)
	fileList = os.listdir(inputDir) 
	for i in fileList:
		os.popen('sox '+inputDir+'/'+i+' '+outputDir+'/'+i+' rate 8k')
		# print('sox '+inputDir+'/'+i+' '+outputDir+'/'+i+' rate 8k')


def Makedir(targetdir):
    if os.path.isdir(targetdir) is False:
        os.makedirs(targetdir)

if __name__ == "__main__":
	main()
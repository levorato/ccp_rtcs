#!/usr/bin/env python
# coding: utf-8

import os
import os.path
import argparse
import subprocess

def main(argv):

    parser = argparse.ArgumentParser(description='Invoke RCCP instance generator.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the result files (one for each experiment)')
    parser.add_argument('--filefilter', default='*_TRACE_*.csv', required=False,
                        help='the file extension for result files (default: *.csv)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    args = parser.parse_args()

    print('Input folders are ', folders)
    print('File filter is ', filter)
    processResult(folders, filter)


#subprocess.check_output(['ls','-l']) #all that is technically needed...
print subprocess.check_output(['ls','-l'])
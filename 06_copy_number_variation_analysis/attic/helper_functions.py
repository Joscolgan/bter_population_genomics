#!/usr/bin/env python3
##############################################################################
# Author: Joe Colgan                   Program: helper_functions.py
#
# Date: 24/07/2015
#
##############################################################################
""" This script contains custom defined functions to:
    1) Check input files are not empty
    2) Check paths for executable tools are true
    3) Check input paired-end files contain the same identifiers"""

# Import modules
import os.path
import tempfile
import subprocess


def check_files_arent_empty(files):
    """ This function takes a list of files as input and for each file present in the list,\
    evaluates if the byte size of the file isn't zero and raises an error if false"""
    for item in files:
        if os.path.getsize(item) == 0:
            raise IOError("File is empty: %s" % (item))


def check_tools(tools):
    """ This function takes a dictionary of paths to software as input and for each file present in\
    the list, evaluates if path is true and raises an error if false"""
    for toolname in tools:
        if not os.path.exists(tools[toolname]):
            raise IOError("Missing tool: %s" % tools[toolname])
        if not os.access(tools[toolname], os.X_OK):
            raise IOError("Cannot execute (try chmod?): %s" % tools[toolname])


def check_identical_fastq_headers(file1, file2):
    """ This function takes two files as input, grabs sequence identifiers and compares to check \
    there is no difference. If there is an difference, an error is raised"""
    mytempfile = tempfile.mkstemp()
    head1 = tempfile.mkstemp()
    head2 = tempfile.mkstemp()

    subprocess.call("grep '@HWI' %s | cut -d ' ' -f 1 - > %s" % (file1, head1[1]), shell=True)
    subprocess.call("grep '@HWI' %s | cut -d ' ' -f 1 - > %s" % (file2, head2[1]), shell=True)

    os.system("diff %s %s > %s" % (head1[1], head2[1], mytempfile[1]))

    if os.path.getsize(mytempfile[1]) != 0:
        raise IOError("Fastq Id differences between: %s, %s" % (file1, file2))

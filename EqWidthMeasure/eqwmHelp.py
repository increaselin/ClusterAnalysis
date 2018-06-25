#! /usr/bin/env python
# 
# Program: eqwmHelp.py
#
# Author: Mike Lum
#
# Usage (Function call): printHelpText()
#
# Description: Help call & instructional text for the eqwm library
#
# Revision History:
#    Date        Vers.    Author        Description
#    2014-06-12    1.0f1    Lum        First checked in

def printHelpText():
    print('Program: eqwm.py\n')
    print('Measures the equivalent widths of the lines passed in the file, as they appear')
    print('in the passed spectra file.')
    print('Output is either placed in <fits spectra file>.log, or in the output file specified')
    print('by the -o flag.')
    print('Usage: ./eqwm.py <fits spectra file> <line list file> [-ahov]\n')
    print('<fits spectra file>: The spectra to measure, in .fits format, regardless of extension.')
    print('<line list file>: A text file, with the following format for each line:')
    print('\tW.L. (air)\tAt.No. + ionization (NN.II)\texcitation pot.\tgf value (fwhm)\tC6 (vanDerWaals) value')
    print('\tThe fields may be separated by any whitespace character(s)')
    print('Options:')
    print('-h: Print this help text.')
    print('-a: Enable automated fitting.')
    print('-c: Include column headers in the output file')
    print('-d: Pause between line measurements when using automated fitting.')
    print('-s: Record unmeasured lines.')
    print('-o <output file name>: Use this file, instead of <fits spectra file>.log')
    print('-v: Use verbose progress and error messages')



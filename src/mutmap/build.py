#!/usr/bin/env python2.7

import sys
import os
import argparse
from subprocess import Popen, PIPE, call

def escapeAsciiForJavaScript(string):
    # Need to escape quotes and newlines for JavaScript
    return string.replace('"', '\\x22').replace('\n', '\\x0A')

def checkRscriptDependencies(rscript):
    for dep in ['kinship2', 'jsonlite']:
        proc = Popen([rscript, '-'], stdout=PIPE, stdin=PIPE,
                     stderr=PIPE)
        proc.communicate(input="library({})".format(dep))

        if proc.returncode != 0:
            print("mutmap requires {} library for R. Please install".format(dep))
            exit(1)

def computeLayout(rscript, sourceDir, pedFileData):

    checkRscriptDependencies(rscript)

    layoutScriptPath = os.path.join(sourceDir, 'pedigree_and_layout.R')
    rProc = Popen([rscript, layoutScriptPath], stdout=PIPE,
                  stdin=PIPE, stderr=PIPE)
    layoutData = rProc.communicate(input=pedFileData)[0]
    return layoutData

def createFileReader(directory):
    def fileReader(filename):
        path = os.path.join(directory, filename)
        with open(path, 'r') as f:
            return f.read()
    return fileReader


if __name__ == '__main__':

    argParser = argparse.ArgumentParser(description="Visualization for denovogear")
    argParser.add_argument("--rscript_location", type=str,
                           required=True,
                           help="Location of Rscript executable")
    argParser.add_argument('-p', "--ped_file", type=str, required=True,
                           help="Pedigree file in ped format")
    argParser.add_argument("--source_dir", type=str, required=True,
                           help="Source directory for input files")
    argParser.add_argument('-i', "--dng_output_file_path", type=str, required=True,
                           help="denovogear output VCF")
    argParser.add_argument('-o', "--output_file_path", type=str, required=True,
                           help="Output file path")
    args = argParser.parse_args()

    readFile = createFileReader(args.source_dir)

    rscript = args.rscript_location

    pedFilePath = args.ped_file

    # Store ped file and compute layout
    with open(pedFilePath, 'r') as pedFile:
        pedFileData = pedFile.read()
    layoutData = computeLayout(rscript, args.source_dir, pedFileData)
    pedFileData = escapeAsciiForJavaScript(pedFileData)

    # Store dng output file
    dngOutputFilePath = args.dng_output_file_path
    with open(dngOutputFilePath, 'r') as dngOutputFile:
        dngOutputData = escapeAsciiForJavaScript(dngOutputFile.read())

    # Load javascript bundle file (from browserify build)
    bundleFileData = readFile('bundle.js')

    # Load template file
    templateFilePath = os.path.join(args.source_dir, 'index.html')
    with open(templateFilePath, 'r') as templateFile:
        templateFileData = templateFile.read()
    outFileData = templateFileData

    # Replace template placeholders with actual data
    outFileData = outFileData.replace(
        '<!--BUNDLE_JS_PLACEHOLDER-->',
        bundleFileData)
    outFileData = outFileData.replace(
        '/*PEDIGREE_FILE_TEXT_PLACEHOLDER*/',
        "var pedigreeFileText = \"" + pedFileData + "\";")
    outFileData = outFileData.replace(
        '/*LAYOUT_DATA_PLACEHOLDER*/',
        "var layoutData = " + layoutData + ";")
    outFileData = outFileData.replace(
        '/*DNG_VCF_DATA_PLACEHOLDER*/',
        "var dngOutputFileText = \"" + dngOutputData + "\";")

    outFileName = args.output_file_path
    with open(outFileName, 'w') as outFile:
        outFile.write(outFileData)

/*
 * Copyright (c) 2010, 2011 Genome Research Ltd.
 * Copyright (c) 2012, 2013 Donald Conrad and Washington University in St. Louis
 * Authors: Donald Conrad <dconrad@genetics.wustl.edu>,
 * Avinash Ramu <aramu@genetics.wustl.edu>
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#define MAX_TRIOS 1000 // max number of trios in the bcf filed
#define MAX_PAIRS 100 // max number of paired samples in the bcf file
#define ID_LENGTH 1000
#define LINE_LENGTH 10000

// Trio Structure
typedef struct  {
  char fID[ID_LENGTH];  // family ID
  char cID[ID_LENGTH];  // child sampleID
  char dID[ID_LENGTH];  // dad sampleID
  char mID[ID_LENGTH];  // mom sampleID
} Trio;

// Pair Structure
typedef struct  {
  char pairID[ID_LENGTH]; // case ID
  char tumorID[ID_LENGTH];  // tumor sample ID
  char normalID[ID_LENGTH]; // normal sample ID
} Pair;

void parse_ped(std::string ped_file, Trio** t, Pair** p, int& trio_count, int& pair_count);

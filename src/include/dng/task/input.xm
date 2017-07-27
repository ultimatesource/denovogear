/*
 * Copyright (c) 2017 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
 *
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

/***************************************************************************
 *    X-Macro List                                                         *
 *                                                                         *
 * XM((long)(name), (shortname), "description", typename, defaultvalue)    *
 ***************************************************************************/

XM((fasta), (f), "faidx indexed reference sequence file", std::string, "")
XM((header), (h), "Location of separate bam header file for read-groups.", 
   std::string, "")
XM((ped), (p), "the pedigree file", std::string, "")
XM((region), (r), "chromosomal region", std::string, "")
XM((rgtag), , "combine read groups using @RG tags, e.g. ID, SM, LB, or DS.",
   std::string, "LB")
XM((sam)(files), (s), "file containing a list of input filenames, one per line",
   std::string, "")

XM((output), (o), "output file", std::string, "-")

XM((min)(qlen), (l), "minimum query length", int, 0)
XM((min)(basequal), (Q), "minimum base quality", int, 13)
XM((min)(mapqual), (q), "minimum mapping quality", int, 0)
XM((normalize)(somatic)(trees), , "scale somatic trees to a height of 1", bool, DL(true, "on"))

/*
 * Copyright (c) 2015 Reed A. Cartwright
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

#include <dng/mutation.h>

const char dng::mitotic_diploid_mutation_labels[10][10][6] = {
    {"AA>AA", "AA>AC", "AA>AG", "AA>AT", "AA>CC", "AA>CG", "AA>CT", "AA>GG", "AA>GT", "AA>TT"},
    {"AC>AA", "AC>AC", "AC>AG", "AC>AT", "AC>CC", "AC>CG", "AC>CT", "AC>GG", "AC>GT", "AC>TT"},
    {"AG>AA", "AG>AC", "AG>AG", "AG>AT", "AG>CC", "AG>CG", "AG>CT", "AG>GG", "AG>GT", "AG>TT"},
    {"AT>AA", "AT>AC", "AT>AG", "AT>AT", "AT>CC", "AT>CG", "AT>CT", "AT>GG", "AT>GT", "AT>TT"},
    {"CC>AA", "CC>AC", "CC>AG", "CC>AT", "CC>CC", "CC>CG", "CC>CT", "CC>GG", "CC>GT", "CC>TT"},
    {"CG>AA", "CG>AC", "CG>AG", "CG>AT", "CG>CC", "CG>CG", "CG>CT", "CG>GG", "CG>GT", "CG>TT"},
    {"CT>AA", "CT>AC", "CT>AG", "CT>AT", "CT>CC", "CT>CG", "CT>CT", "CT>GG", "CT>GT", "CT>TT"},
    {"GG>AA", "GG>AC", "GG>AG", "GG>AT", "GG>CC", "GG>CG", "GG>CT", "GG>GG", "GG>GT", "GG>TT"},
    {"GT>AA", "GT>AC", "GT>AG", "GT>AT", "GT>CC", "GT>CG", "GT>CT", "GT>GG", "GT>GT", "GT>TT"},
    {"TT>AA", "TT>AC", "TT>AG", "TT>AT", "TT>CC", "TT>CG", "TT>CT", "TT>GG", "TT>GT", "TT>TT"}
};

const char dng::meiotic_diploid_mutation_labels[100][10][9] = {
    {"AAxAA>AA", "AAxAA>AC", "AAxAA>AG", "AAxAA>AT", "AAxAA>CC", "AAxAA>CG", "AAxAA>CT", "AAxAA>GG", "AAxAA>GT", "AAxAA>TT"},
    {"AAxAC>AA", "AAxAC>AC", "AAxAC>AG", "AAxAC>AT", "AAxAC>CC", "AAxAC>CG", "AAxAC>CT", "AAxAC>GG", "AAxAC>GT", "AAxAC>TT"},
    {"AAxAG>AA", "AAxAG>AC", "AAxAG>AG", "AAxAG>AT", "AAxAG>CC", "AAxAG>CG", "AAxAG>CT", "AAxAG>GG", "AAxAG>GT", "AAxAG>TT"},
    {"AAxAT>AA", "AAxAT>AC", "AAxAT>AG", "AAxAT>AT", "AAxAT>CC", "AAxAT>CG", "AAxAT>CT", "AAxAT>GG", "AAxAT>GT", "AAxAT>TT"},
    {"AAxCC>AA", "AAxCC>AC", "AAxCC>AG", "AAxCC>AT", "AAxCC>CC", "AAxCC>CG", "AAxCC>CT", "AAxCC>GG", "AAxCC>GT", "AAxCC>TT"},
    {"AAxCG>AA", "AAxCG>AC", "AAxCG>AG", "AAxCG>AT", "AAxCG>CC", "AAxCG>CG", "AAxCG>CT", "AAxCG>GG", "AAxCG>GT", "AAxCG>TT"},
    {"AAxCT>AA", "AAxCT>AC", "AAxCT>AG", "AAxCT>AT", "AAxCT>CC", "AAxCT>CG", "AAxCT>CT", "AAxCT>GG", "AAxCT>GT", "AAxCT>TT"},
    {"AAxGG>AA", "AAxGG>AC", "AAxGG>AG", "AAxGG>AT", "AAxGG>CC", "AAxGG>CG", "AAxGG>CT", "AAxGG>GG", "AAxGG>GT", "AAxGG>TT"},
    {"AAxGT>AA", "AAxGT>AC", "AAxGT>AG", "AAxGT>AT", "AAxGT>CC", "AAxGT>CG", "AAxGT>CT", "AAxGT>GG", "AAxGT>GT", "AAxGT>TT"},
    {"AAxTT>AA", "AAxTT>AC", "AAxTT>AG", "AAxTT>AT", "AAxTT>CC", "AAxTT>CG", "AAxTT>CT", "AAxTT>GG", "AAxTT>GT", "AAxTT>TT"},
    {"ACxAA>AA", "ACxAA>AC", "ACxAA>AG", "ACxAA>AT", "ACxAA>CC", "ACxAA>CG", "ACxAA>CT", "ACxAA>GG", "ACxAA>GT", "ACxAA>TT"},
    {"ACxAC>AA", "ACxAC>AC", "ACxAC>AG", "ACxAC>AT", "ACxAC>CC", "ACxAC>CG", "ACxAC>CT", "ACxAC>GG", "ACxAC>GT", "ACxAC>TT"},
    {"ACxAG>AA", "ACxAG>AC", "ACxAG>AG", "ACxAG>AT", "ACxAG>CC", "ACxAG>CG", "ACxAG>CT", "ACxAG>GG", "ACxAG>GT", "ACxAG>TT"},
    {"ACxAT>AA", "ACxAT>AC", "ACxAT>AG", "ACxAT>AT", "ACxAT>CC", "ACxAT>CG", "ACxAT>CT", "ACxAT>GG", "ACxAT>GT", "ACxAT>TT"},
    {"ACxCC>AA", "ACxCC>AC", "ACxCC>AG", "ACxCC>AT", "ACxCC>CC", "ACxCC>CG", "ACxCC>CT", "ACxCC>GG", "ACxCC>GT", "ACxCC>TT"},
    {"ACxCG>AA", "ACxCG>AC", "ACxCG>AG", "ACxCG>AT", "ACxCG>CC", "ACxCG>CG", "ACxCG>CT", "ACxCG>GG", "ACxCG>GT", "ACxCG>TT"},
    {"ACxCT>AA", "ACxCT>AC", "ACxCT>AG", "ACxCT>AT", "ACxCT>CC", "ACxCT>CG", "ACxCT>CT", "ACxCT>GG", "ACxCT>GT", "ACxCT>TT"},
    {"ACxGG>AA", "ACxGG>AC", "ACxGG>AG", "ACxGG>AT", "ACxGG>CC", "ACxGG>CG", "ACxGG>CT", "ACxGG>GG", "ACxGG>GT", "ACxGG>TT"},
    {"ACxGT>AA", "ACxGT>AC", "ACxGT>AG", "ACxGT>AT", "ACxGT>CC", "ACxGT>CG", "ACxGT>CT", "ACxGT>GG", "ACxGT>GT", "ACxGT>TT"},
    {"ACxTT>AA", "ACxTT>AC", "ACxTT>AG", "ACxTT>AT", "ACxTT>CC", "ACxTT>CG", "ACxTT>CT", "ACxTT>GG", "ACxTT>GT", "ACxTT>TT"},
    {"AGxAA>AA", "AGxAA>AC", "AGxAA>AG", "AGxAA>AT", "AGxAA>CC", "AGxAA>CG", "AGxAA>CT", "AGxAA>GG", "AGxAA>GT", "AGxAA>TT"},
    {"AGxAC>AA", "AGxAC>AC", "AGxAC>AG", "AGxAC>AT", "AGxAC>CC", "AGxAC>CG", "AGxAC>CT", "AGxAC>GG", "AGxAC>GT", "AGxAC>TT"},
    {"AGxAG>AA", "AGxAG>AC", "AGxAG>AG", "AGxAG>AT", "AGxAG>CC", "AGxAG>CG", "AGxAG>CT", "AGxAG>GG", "AGxAG>GT", "AGxAG>TT"},
    {"AGxAT>AA", "AGxAT>AC", "AGxAT>AG", "AGxAT>AT", "AGxAT>CC", "AGxAT>CG", "AGxAT>CT", "AGxAT>GG", "AGxAT>GT", "AGxAT>TT"},
    {"AGxCC>AA", "AGxCC>AC", "AGxCC>AG", "AGxCC>AT", "AGxCC>CC", "AGxCC>CG", "AGxCC>CT", "AGxCC>GG", "AGxCC>GT", "AGxCC>TT"},
    {"AGxCG>AA", "AGxCG>AC", "AGxCG>AG", "AGxCG>AT", "AGxCG>CC", "AGxCG>CG", "AGxCG>CT", "AGxCG>GG", "AGxCG>GT", "AGxCG>TT"},
    {"AGxCT>AA", "AGxCT>AC", "AGxCT>AG", "AGxCT>AT", "AGxCT>CC", "AGxCT>CG", "AGxCT>CT", "AGxCT>GG", "AGxCT>GT", "AGxCT>TT"},
    {"AGxGG>AA", "AGxGG>AC", "AGxGG>AG", "AGxGG>AT", "AGxGG>CC", "AGxGG>CG", "AGxGG>CT", "AGxGG>GG", "AGxGG>GT", "AGxGG>TT"},
    {"AGxGT>AA", "AGxGT>AC", "AGxGT>AG", "AGxGT>AT", "AGxGT>CC", "AGxGT>CG", "AGxGT>CT", "AGxGT>GG", "AGxGT>GT", "AGxGT>TT"},
    {"AGxTT>AA", "AGxTT>AC", "AGxTT>AG", "AGxTT>AT", "AGxTT>CC", "AGxTT>CG", "AGxTT>CT", "AGxTT>GG", "AGxTT>GT", "AGxTT>TT"},
    {"ATxAA>AA", "ATxAA>AC", "ATxAA>AG", "ATxAA>AT", "ATxAA>CC", "ATxAA>CG", "ATxAA>CT", "ATxAA>GG", "ATxAA>GT", "ATxAA>TT"},
    {"ATxAC>AA", "ATxAC>AC", "ATxAC>AG", "ATxAC>AT", "ATxAC>CC", "ATxAC>CG", "ATxAC>CT", "ATxAC>GG", "ATxAC>GT", "ATxAC>TT"},
    {"ATxAG>AA", "ATxAG>AC", "ATxAG>AG", "ATxAG>AT", "ATxAG>CC", "ATxAG>CG", "ATxAG>CT", "ATxAG>GG", "ATxAG>GT", "ATxAG>TT"},
    {"ATxAT>AA", "ATxAT>AC", "ATxAT>AG", "ATxAT>AT", "ATxAT>CC", "ATxAT>CG", "ATxAT>CT", "ATxAT>GG", "ATxAT>GT", "ATxAT>TT"},
    {"ATxCC>AA", "ATxCC>AC", "ATxCC>AG", "ATxCC>AT", "ATxCC>CC", "ATxCC>CG", "ATxCC>CT", "ATxCC>GG", "ATxCC>GT", "ATxCC>TT"},
    {"ATxCG>AA", "ATxCG>AC", "ATxCG>AG", "ATxCG>AT", "ATxCG>CC", "ATxCG>CG", "ATxCG>CT", "ATxCG>GG", "ATxCG>GT", "ATxCG>TT"},
    {"ATxCT>AA", "ATxCT>AC", "ATxCT>AG", "ATxCT>AT", "ATxCT>CC", "ATxCT>CG", "ATxCT>CT", "ATxCT>GG", "ATxCT>GT", "ATxCT>TT"},
    {"ATxGG>AA", "ATxGG>AC", "ATxGG>AG", "ATxGG>AT", "ATxGG>CC", "ATxGG>CG", "ATxGG>CT", "ATxGG>GG", "ATxGG>GT", "ATxGG>TT"},
    {"ATxGT>AA", "ATxGT>AC", "ATxGT>AG", "ATxGT>AT", "ATxGT>CC", "ATxGT>CG", "ATxGT>CT", "ATxGT>GG", "ATxGT>GT", "ATxGT>TT"},
    {"ATxTT>AA", "ATxTT>AC", "ATxTT>AG", "ATxTT>AT", "ATxTT>CC", "ATxTT>CG", "ATxTT>CT", "ATxTT>GG", "ATxTT>GT", "ATxTT>TT"},
    {"CCxAA>AA", "CCxAA>AC", "CCxAA>AG", "CCxAA>AT", "CCxAA>CC", "CCxAA>CG", "CCxAA>CT", "CCxAA>GG", "CCxAA>GT", "CCxAA>TT"},
    {"CCxAC>AA", "CCxAC>AC", "CCxAC>AG", "CCxAC>AT", "CCxAC>CC", "CCxAC>CG", "CCxAC>CT", "CCxAC>GG", "CCxAC>GT", "CCxAC>TT"},
    {"CCxAG>AA", "CCxAG>AC", "CCxAG>AG", "CCxAG>AT", "CCxAG>CC", "CCxAG>CG", "CCxAG>CT", "CCxAG>GG", "CCxAG>GT", "CCxAG>TT"},
    {"CCxAT>AA", "CCxAT>AC", "CCxAT>AG", "CCxAT>AT", "CCxAT>CC", "CCxAT>CG", "CCxAT>CT", "CCxAT>GG", "CCxAT>GT", "CCxAT>TT"},
    {"CCxCC>AA", "CCxCC>AC", "CCxCC>AG", "CCxCC>AT", "CCxCC>CC", "CCxCC>CG", "CCxCC>CT", "CCxCC>GG", "CCxCC>GT", "CCxCC>TT"},
    {"CCxCG>AA", "CCxCG>AC", "CCxCG>AG", "CCxCG>AT", "CCxCG>CC", "CCxCG>CG", "CCxCG>CT", "CCxCG>GG", "CCxCG>GT", "CCxCG>TT"},
    {"CCxCT>AA", "CCxCT>AC", "CCxCT>AG", "CCxCT>AT", "CCxCT>CC", "CCxCT>CG", "CCxCT>CT", "CCxCT>GG", "CCxCT>GT", "CCxCT>TT"},
    {"CCxGG>AA", "CCxGG>AC", "CCxGG>AG", "CCxGG>AT", "CCxGG>CC", "CCxGG>CG", "CCxGG>CT", "CCxGG>GG", "CCxGG>GT", "CCxGG>TT"},
    {"CCxGT>AA", "CCxGT>AC", "CCxGT>AG", "CCxGT>AT", "CCxGT>CC", "CCxGT>CG", "CCxGT>CT", "CCxGT>GG", "CCxGT>GT", "CCxGT>TT"},
    {"CCxTT>AA", "CCxTT>AC", "CCxTT>AG", "CCxTT>AT", "CCxTT>CC", "CCxTT>CG", "CCxTT>CT", "CCxTT>GG", "CCxTT>GT", "CCxTT>TT"},
    {"CGxAA>AA", "CGxAA>AC", "CGxAA>AG", "CGxAA>AT", "CGxAA>CC", "CGxAA>CG", "CGxAA>CT", "CGxAA>GG", "CGxAA>GT", "CGxAA>TT"},
    {"CGxAC>AA", "CGxAC>AC", "CGxAC>AG", "CGxAC>AT", "CGxAC>CC", "CGxAC>CG", "CGxAC>CT", "CGxAC>GG", "CGxAC>GT", "CGxAC>TT"},
    {"CGxAG>AA", "CGxAG>AC", "CGxAG>AG", "CGxAG>AT", "CGxAG>CC", "CGxAG>CG", "CGxAG>CT", "CGxAG>GG", "CGxAG>GT", "CGxAG>TT"},
    {"CGxAT>AA", "CGxAT>AC", "CGxAT>AG", "CGxAT>AT", "CGxAT>CC", "CGxAT>CG", "CGxAT>CT", "CGxAT>GG", "CGxAT>GT", "CGxAT>TT"},
    {"CGxCC>AA", "CGxCC>AC", "CGxCC>AG", "CGxCC>AT", "CGxCC>CC", "CGxCC>CG", "CGxCC>CT", "CGxCC>GG", "CGxCC>GT", "CGxCC>TT"},
    {"CGxCG>AA", "CGxCG>AC", "CGxCG>AG", "CGxCG>AT", "CGxCG>CC", "CGxCG>CG", "CGxCG>CT", "CGxCG>GG", "CGxCG>GT", "CGxCG>TT"},
    {"CGxCT>AA", "CGxCT>AC", "CGxCT>AG", "CGxCT>AT", "CGxCT>CC", "CGxCT>CG", "CGxCT>CT", "CGxCT>GG", "CGxCT>GT", "CGxCT>TT"},
    {"CGxGG>AA", "CGxGG>AC", "CGxGG>AG", "CGxGG>AT", "CGxGG>CC", "CGxGG>CG", "CGxGG>CT", "CGxGG>GG", "CGxGG>GT", "CGxGG>TT"},
    {"CGxGT>AA", "CGxGT>AC", "CGxGT>AG", "CGxGT>AT", "CGxGT>CC", "CGxGT>CG", "CGxGT>CT", "CGxGT>GG", "CGxGT>GT", "CGxGT>TT"},
    {"CGxTT>AA", "CGxTT>AC", "CGxTT>AG", "CGxTT>AT", "CGxTT>CC", "CGxTT>CG", "CGxTT>CT", "CGxTT>GG", "CGxTT>GT", "CGxTT>TT"},
    {"CTxAA>AA", "CTxAA>AC", "CTxAA>AG", "CTxAA>AT", "CTxAA>CC", "CTxAA>CG", "CTxAA>CT", "CTxAA>GG", "CTxAA>GT", "CTxAA>TT"},
    {"CTxAC>AA", "CTxAC>AC", "CTxAC>AG", "CTxAC>AT", "CTxAC>CC", "CTxAC>CG", "CTxAC>CT", "CTxAC>GG", "CTxAC>GT", "CTxAC>TT"},
    {"CTxAG>AA", "CTxAG>AC", "CTxAG>AG", "CTxAG>AT", "CTxAG>CC", "CTxAG>CG", "CTxAG>CT", "CTxAG>GG", "CTxAG>GT", "CTxAG>TT"},
    {"CTxAT>AA", "CTxAT>AC", "CTxAT>AG", "CTxAT>AT", "CTxAT>CC", "CTxAT>CG", "CTxAT>CT", "CTxAT>GG", "CTxAT>GT", "CTxAT>TT"},
    {"CTxCC>AA", "CTxCC>AC", "CTxCC>AG", "CTxCC>AT", "CTxCC>CC", "CTxCC>CG", "CTxCC>CT", "CTxCC>GG", "CTxCC>GT", "CTxCC>TT"},
    {"CTxCG>AA", "CTxCG>AC", "CTxCG>AG", "CTxCG>AT", "CTxCG>CC", "CTxCG>CG", "CTxCG>CT", "CTxCG>GG", "CTxCG>GT", "CTxCG>TT"},
    {"CTxCT>AA", "CTxCT>AC", "CTxCT>AG", "CTxCT>AT", "CTxCT>CC", "CTxCT>CG", "CTxCT>CT", "CTxCT>GG", "CTxCT>GT", "CTxCT>TT"},
    {"CTxGG>AA", "CTxGG>AC", "CTxGG>AG", "CTxGG>AT", "CTxGG>CC", "CTxGG>CG", "CTxGG>CT", "CTxGG>GG", "CTxGG>GT", "CTxGG>TT"},
    {"CTxGT>AA", "CTxGT>AC", "CTxGT>AG", "CTxGT>AT", "CTxGT>CC", "CTxGT>CG", "CTxGT>CT", "CTxGT>GG", "CTxGT>GT", "CTxGT>TT"},
    {"CTxTT>AA", "CTxTT>AC", "CTxTT>AG", "CTxTT>AT", "CTxTT>CC", "CTxTT>CG", "CTxTT>CT", "CTxTT>GG", "CTxTT>GT", "CTxTT>TT"},
    {"GGxAA>AA", "GGxAA>AC", "GGxAA>AG", "GGxAA>AT", "GGxAA>CC", "GGxAA>CG", "GGxAA>CT", "GGxAA>GG", "GGxAA>GT", "GGxAA>TT"},
    {"GGxAC>AA", "GGxAC>AC", "GGxAC>AG", "GGxAC>AT", "GGxAC>CC", "GGxAC>CG", "GGxAC>CT", "GGxAC>GG", "GGxAC>GT", "GGxAC>TT"},
    {"GGxAG>AA", "GGxAG>AC", "GGxAG>AG", "GGxAG>AT", "GGxAG>CC", "GGxAG>CG", "GGxAG>CT", "GGxAG>GG", "GGxAG>GT", "GGxAG>TT"},
    {"GGxAT>AA", "GGxAT>AC", "GGxAT>AG", "GGxAT>AT", "GGxAT>CC", "GGxAT>CG", "GGxAT>CT", "GGxAT>GG", "GGxAT>GT", "GGxAT>TT"},
    {"GGxCC>AA", "GGxCC>AC", "GGxCC>AG", "GGxCC>AT", "GGxCC>CC", "GGxCC>CG", "GGxCC>CT", "GGxCC>GG", "GGxCC>GT", "GGxCC>TT"},
    {"GGxCG>AA", "GGxCG>AC", "GGxCG>AG", "GGxCG>AT", "GGxCG>CC", "GGxCG>CG", "GGxCG>CT", "GGxCG>GG", "GGxCG>GT", "GGxCG>TT"},
    {"GGxCT>AA", "GGxCT>AC", "GGxCT>AG", "GGxCT>AT", "GGxCT>CC", "GGxCT>CG", "GGxCT>CT", "GGxCT>GG", "GGxCT>GT", "GGxCT>TT"},
    {"GGxGG>AA", "GGxGG>AC", "GGxGG>AG", "GGxGG>AT", "GGxGG>CC", "GGxGG>CG", "GGxGG>CT", "GGxGG>GG", "GGxGG>GT", "GGxGG>TT"},
    {"GGxGT>AA", "GGxGT>AC", "GGxGT>AG", "GGxGT>AT", "GGxGT>CC", "GGxGT>CG", "GGxGT>CT", "GGxGT>GG", "GGxGT>GT", "GGxGT>TT"},
    {"GGxTT>AA", "GGxTT>AC", "GGxTT>AG", "GGxTT>AT", "GGxTT>CC", "GGxTT>CG", "GGxTT>CT", "GGxTT>GG", "GGxTT>GT", "GGxTT>TT"},
    {"GTxAA>AA", "GTxAA>AC", "GTxAA>AG", "GTxAA>AT", "GTxAA>CC", "GTxAA>CG", "GTxAA>CT", "GTxAA>GG", "GTxAA>GT", "GTxAA>TT"},
    {"GTxAC>AA", "GTxAC>AC", "GTxAC>AG", "GTxAC>AT", "GTxAC>CC", "GTxAC>CG", "GTxAC>CT", "GTxAC>GG", "GTxAC>GT", "GTxAC>TT"},
    {"GTxAG>AA", "GTxAG>AC", "GTxAG>AG", "GTxAG>AT", "GTxAG>CC", "GTxAG>CG", "GTxAG>CT", "GTxAG>GG", "GTxAG>GT", "GTxAG>TT"},
    {"GTxAT>AA", "GTxAT>AC", "GTxAT>AG", "GTxAT>AT", "GTxAT>CC", "GTxAT>CG", "GTxAT>CT", "GTxAT>GG", "GTxAT>GT", "GTxAT>TT"},
    {"GTxCC>AA", "GTxCC>AC", "GTxCC>AG", "GTxCC>AT", "GTxCC>CC", "GTxCC>CG", "GTxCC>CT", "GTxCC>GG", "GTxCC>GT", "GTxCC>TT"},
    {"GTxCG>AA", "GTxCG>AC", "GTxCG>AG", "GTxCG>AT", "GTxCG>CC", "GTxCG>CG", "GTxCG>CT", "GTxCG>GG", "GTxCG>GT", "GTxCG>TT"},
    {"GTxCT>AA", "GTxCT>AC", "GTxCT>AG", "GTxCT>AT", "GTxCT>CC", "GTxCT>CG", "GTxCT>CT", "GTxCT>GG", "GTxCT>GT", "GTxCT>TT"},
    {"GTxGG>AA", "GTxGG>AC", "GTxGG>AG", "GTxGG>AT", "GTxGG>CC", "GTxGG>CG", "GTxGG>CT", "GTxGG>GG", "GTxGG>GT", "GTxGG>TT"},
    {"GTxGT>AA", "GTxGT>AC", "GTxGT>AG", "GTxGT>AT", "GTxGT>CC", "GTxGT>CG", "GTxGT>CT", "GTxGT>GG", "GTxGT>GT", "GTxGT>TT"},
    {"GTxTT>AA", "GTxTT>AC", "GTxTT>AG", "GTxTT>AT", "GTxTT>CC", "GTxTT>CG", "GTxTT>CT", "GTxTT>GG", "GTxTT>GT", "GTxTT>TT"},
    {"TTxAA>AA", "TTxAA>AC", "TTxAA>AG", "TTxAA>AT", "TTxAA>CC", "TTxAA>CG", "TTxAA>CT", "TTxAA>GG", "TTxAA>GT", "TTxAA>TT"},
    {"TTxAC>AA", "TTxAC>AC", "TTxAC>AG", "TTxAC>AT", "TTxAC>CC", "TTxAC>CG", "TTxAC>CT", "TTxAC>GG", "TTxAC>GT", "TTxAC>TT"},
    {"TTxAG>AA", "TTxAG>AC", "TTxAG>AG", "TTxAG>AT", "TTxAG>CC", "TTxAG>CG", "TTxAG>CT", "TTxAG>GG", "TTxAG>GT", "TTxAG>TT"},
    {"TTxAT>AA", "TTxAT>AC", "TTxAT>AG", "TTxAT>AT", "TTxAT>CC", "TTxAT>CG", "TTxAT>CT", "TTxAT>GG", "TTxAT>GT", "TTxAT>TT"},
    {"TTxCC>AA", "TTxCC>AC", "TTxCC>AG", "TTxCC>AT", "TTxCC>CC", "TTxCC>CG", "TTxCC>CT", "TTxCC>GG", "TTxCC>GT", "TTxCC>TT"},
    {"TTxCG>AA", "TTxCG>AC", "TTxCG>AG", "TTxCG>AT", "TTxCG>CC", "TTxCG>CG", "TTxCG>CT", "TTxCG>GG", "TTxCG>GT", "TTxCG>TT"},
    {"TTxCT>AA", "TTxCT>AC", "TTxCT>AG", "TTxCT>AT", "TTxCT>CC", "TTxCT>CG", "TTxCT>CT", "TTxCT>GG", "TTxCT>GT", "TTxCT>TT"},
    {"TTxGG>AA", "TTxGG>AC", "TTxGG>AG", "TTxGG>AT", "TTxGG>CC", "TTxGG>CG", "TTxGG>CT", "TTxGG>GG", "TTxGG>GT", "TTxGG>TT"},
    {"TTxGT>AA", "TTxGT>AC", "TTxGT>AG", "TTxGT>AT", "TTxGT>CC", "TTxGT>CG", "TTxGT>CT", "TTxGT>GG", "TTxGT>GT", "TTxGT>TT"},
    {"TTxTT>AA", "TTxTT>AC", "TTxTT>AG", "TTxTT>AT", "TTxTT>CC", "TTxTT>CG", "TTxTT>CT", "TTxTT>GG", "TTxTT>GT", "TTxTT>TT"}
};
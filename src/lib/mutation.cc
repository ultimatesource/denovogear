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
    {"AA>AA", "AA>CC", "AA>GG", "AA>TT", "AA>AC", "AA>AG", "AA>AT", "AA>CG", "AA>CT", "AA>GT"},
    {"CC>AA", "CC>CC", "CC>GG", "CC>TT", "CC>AC", "CC>AG", "CC>AT", "CC>CG", "CC>CT", "CC>GT"},
    {"GG>AA", "GG>CC", "GG>GG", "GG>TT", "GG>AC", "GG>AG", "GG>AT", "GG>CG", "GG>CT", "GG>GT"},
    {"TT>AA", "TT>CC", "TT>GG", "TT>TT", "TT>AC", "TT>AG", "TT>AT", "TT>CG", "TT>CT", "TT>GT"},
    {"AC>AA", "AC>CC", "AC>GG", "AC>TT", "AC>AC", "AC>AG", "AC>AT", "AC>CG", "AC>CT", "AC>GT"},
    {"AG>AA", "AG>CC", "AG>GG", "AG>TT", "AG>AC", "AG>AG", "AG>AT", "AG>CG", "AG>CT", "AG>GT"},
    {"AT>AA", "AT>CC", "AT>GG", "AT>TT", "AT>AC", "AT>AG", "AT>AT", "AT>CG", "AT>CT", "AT>GT"},
    {"CG>AA", "CG>CC", "CG>GG", "CG>TT", "CG>AC", "CG>AG", "CG>AT", "CG>CG", "CG>CT", "CG>GT"},
    {"CT>AA", "CT>CC", "CT>GG", "CT>TT", "CT>AC", "CT>AG", "CT>AT", "CT>CG", "CT>CT", "CT>GT"},
    {"GT>AA", "GT>CC", "GT>GG", "GT>TT", "GT>AC", "GT>AG", "GT>AT", "GT>CG", "GT>CT", "GT>GT"}
};
const char dng::meiotic_diploid_mutation_labels[100][10][9] = {
    {"AAxAA>AA", "AAxAA>CC", "AAxAA>GG", "AAxAA>TT", "AAxAA>AC", "AAxAA>AG", "AAxAA>AT", "AAxAA>CG", "AAxAA>CT", "AAxAA>GT"},
    {"AAxCC>AA", "AAxCC>CC", "AAxCC>GG", "AAxCC>TT", "AAxCC>AC", "AAxCC>AG", "AAxCC>AT", "AAxCC>CG", "AAxCC>CT", "AAxCC>GT"},
    {"AAxGG>AA", "AAxGG>CC", "AAxGG>GG", "AAxGG>TT", "AAxGG>AC", "AAxGG>AG", "AAxGG>AT", "AAxGG>CG", "AAxGG>CT", "AAxGG>GT"},
    {"AAxTT>AA", "AAxTT>CC", "AAxTT>GG", "AAxTT>TT", "AAxTT>AC", "AAxTT>AG", "AAxTT>AT", "AAxTT>CG", "AAxTT>CT", "AAxTT>GT"},
    {"AAxAC>AA", "AAxAC>CC", "AAxAC>GG", "AAxAC>TT", "AAxAC>AC", "AAxAC>AG", "AAxAC>AT", "AAxAC>CG", "AAxAC>CT", "AAxAC>GT"},
    {"AAxAG>AA", "AAxAG>CC", "AAxAG>GG", "AAxAG>TT", "AAxAG>AC", "AAxAG>AG", "AAxAG>AT", "AAxAG>CG", "AAxAG>CT", "AAxAG>GT"},
    {"AAxAT>AA", "AAxAT>CC", "AAxAT>GG", "AAxAT>TT", "AAxAT>AC", "AAxAT>AG", "AAxAT>AT", "AAxAT>CG", "AAxAT>CT", "AAxAT>GT"},
    {"AAxCG>AA", "AAxCG>CC", "AAxCG>GG", "AAxCG>TT", "AAxCG>AC", "AAxCG>AG", "AAxCG>AT", "AAxCG>CG", "AAxCG>CT", "AAxCG>GT"},
    {"AAxCT>AA", "AAxCT>CC", "AAxCT>GG", "AAxCT>TT", "AAxCT>AC", "AAxCT>AG", "AAxCT>AT", "AAxCT>CG", "AAxCT>CT", "AAxCT>GT"},
    {"AAxGT>AA", "AAxGT>CC", "AAxGT>GG", "AAxGT>TT", "AAxGT>AC", "AAxGT>AG", "AAxGT>AT", "AAxGT>CG", "AAxGT>CT", "AAxGT>GT"},
    {"CCxAA>AA", "CCxAA>CC", "CCxAA>GG", "CCxAA>TT", "CCxAA>AC", "CCxAA>AG", "CCxAA>AT", "CCxAA>CG", "CCxAA>CT", "CCxAA>GT"},
    {"CCxCC>AA", "CCxCC>CC", "CCxCC>GG", "CCxCC>TT", "CCxCC>AC", "CCxCC>AG", "CCxCC>AT", "CCxCC>CG", "CCxCC>CT", "CCxCC>GT"},
    {"CCxGG>AA", "CCxGG>CC", "CCxGG>GG", "CCxGG>TT", "CCxGG>AC", "CCxGG>AG", "CCxGG>AT", "CCxGG>CG", "CCxGG>CT", "CCxGG>GT"},
    {"CCxTT>AA", "CCxTT>CC", "CCxTT>GG", "CCxTT>TT", "CCxTT>AC", "CCxTT>AG", "CCxTT>AT", "CCxTT>CG", "CCxTT>CT", "CCxTT>GT"},
    {"CCxAC>AA", "CCxAC>CC", "CCxAC>GG", "CCxAC>TT", "CCxAC>AC", "CCxAC>AG", "CCxAC>AT", "CCxAC>CG", "CCxAC>CT", "CCxAC>GT"},
    {"CCxAG>AA", "CCxAG>CC", "CCxAG>GG", "CCxAG>TT", "CCxAG>AC", "CCxAG>AG", "CCxAG>AT", "CCxAG>CG", "CCxAG>CT", "CCxAG>GT"},
    {"CCxAT>AA", "CCxAT>CC", "CCxAT>GG", "CCxAT>TT", "CCxAT>AC", "CCxAT>AG", "CCxAT>AT", "CCxAT>CG", "CCxAT>CT", "CCxAT>GT"},
    {"CCxCG>AA", "CCxCG>CC", "CCxCG>GG", "CCxCG>TT", "CCxCG>AC", "CCxCG>AG", "CCxCG>AT", "CCxCG>CG", "CCxCG>CT", "CCxCG>GT"},
    {"CCxCT>AA", "CCxCT>CC", "CCxCT>GG", "CCxCT>TT", "CCxCT>AC", "CCxCT>AG", "CCxCT>AT", "CCxCT>CG", "CCxCT>CT", "CCxCT>GT"},
    {"CCxGT>AA", "CCxGT>CC", "CCxGT>GG", "CCxGT>TT", "CCxGT>AC", "CCxGT>AG", "CCxGT>AT", "CCxGT>CG", "CCxGT>CT", "CCxGT>GT"},
    {"GGxAA>AA", "GGxAA>CC", "GGxAA>GG", "GGxAA>TT", "GGxAA>AC", "GGxAA>AG", "GGxAA>AT", "GGxAA>CG", "GGxAA>CT", "GGxAA>GT"},
    {"GGxCC>AA", "GGxCC>CC", "GGxCC>GG", "GGxCC>TT", "GGxCC>AC", "GGxCC>AG", "GGxCC>AT", "GGxCC>CG", "GGxCC>CT", "GGxCC>GT"},
    {"GGxGG>AA", "GGxGG>CC", "GGxGG>GG", "GGxGG>TT", "GGxGG>AC", "GGxGG>AG", "GGxGG>AT", "GGxGG>CG", "GGxGG>CT", "GGxGG>GT"},
    {"GGxTT>AA", "GGxTT>CC", "GGxTT>GG", "GGxTT>TT", "GGxTT>AC", "GGxTT>AG", "GGxTT>AT", "GGxTT>CG", "GGxTT>CT", "GGxTT>GT"},
    {"GGxAC>AA", "GGxAC>CC", "GGxAC>GG", "GGxAC>TT", "GGxAC>AC", "GGxAC>AG", "GGxAC>AT", "GGxAC>CG", "GGxAC>CT", "GGxAC>GT"},
    {"GGxAG>AA", "GGxAG>CC", "GGxAG>GG", "GGxAG>TT", "GGxAG>AC", "GGxAG>AG", "GGxAG>AT", "GGxAG>CG", "GGxAG>CT", "GGxAG>GT"},
    {"GGxAT>AA", "GGxAT>CC", "GGxAT>GG", "GGxAT>TT", "GGxAT>AC", "GGxAT>AG", "GGxAT>AT", "GGxAT>CG", "GGxAT>CT", "GGxAT>GT"},
    {"GGxCG>AA", "GGxCG>CC", "GGxCG>GG", "GGxCG>TT", "GGxCG>AC", "GGxCG>AG", "GGxCG>AT", "GGxCG>CG", "GGxCG>CT", "GGxCG>GT"},
    {"GGxCT>AA", "GGxCT>CC", "GGxCT>GG", "GGxCT>TT", "GGxCT>AC", "GGxCT>AG", "GGxCT>AT", "GGxCT>CG", "GGxCT>CT", "GGxCT>GT"},
    {"GGxGT>AA", "GGxGT>CC", "GGxGT>GG", "GGxGT>TT", "GGxGT>AC", "GGxGT>AG", "GGxGT>AT", "GGxGT>CG", "GGxGT>CT", "GGxGT>GT"},
    {"TTxAA>AA", "TTxAA>CC", "TTxAA>GG", "TTxAA>TT", "TTxAA>AC", "TTxAA>AG", "TTxAA>AT", "TTxAA>CG", "TTxAA>CT", "TTxAA>GT"},
    {"TTxCC>AA", "TTxCC>CC", "TTxCC>GG", "TTxCC>TT", "TTxCC>AC", "TTxCC>AG", "TTxCC>AT", "TTxCC>CG", "TTxCC>CT", "TTxCC>GT"},
    {"TTxGG>AA", "TTxGG>CC", "TTxGG>GG", "TTxGG>TT", "TTxGG>AC", "TTxGG>AG", "TTxGG>AT", "TTxGG>CG", "TTxGG>CT", "TTxGG>GT"},
    {"TTxTT>AA", "TTxTT>CC", "TTxTT>GG", "TTxTT>TT", "TTxTT>AC", "TTxTT>AG", "TTxTT>AT", "TTxTT>CG", "TTxTT>CT", "TTxTT>GT"},
    {"TTxAC>AA", "TTxAC>CC", "TTxAC>GG", "TTxAC>TT", "TTxAC>AC", "TTxAC>AG", "TTxAC>AT", "TTxAC>CG", "TTxAC>CT", "TTxAC>GT"},
    {"TTxAG>AA", "TTxAG>CC", "TTxAG>GG", "TTxAG>TT", "TTxAG>AC", "TTxAG>AG", "TTxAG>AT", "TTxAG>CG", "TTxAG>CT", "TTxAG>GT"},
    {"TTxAT>AA", "TTxAT>CC", "TTxAT>GG", "TTxAT>TT", "TTxAT>AC", "TTxAT>AG", "TTxAT>AT", "TTxAT>CG", "TTxAT>CT", "TTxAT>GT"},
    {"TTxCG>AA", "TTxCG>CC", "TTxCG>GG", "TTxCG>TT", "TTxCG>AC", "TTxCG>AG", "TTxCG>AT", "TTxCG>CG", "TTxCG>CT", "TTxCG>GT"},
    {"TTxCT>AA", "TTxCT>CC", "TTxCT>GG", "TTxCT>TT", "TTxCT>AC", "TTxCT>AG", "TTxCT>AT", "TTxCT>CG", "TTxCT>CT", "TTxCT>GT"},
    {"TTxGT>AA", "TTxGT>CC", "TTxGT>GG", "TTxGT>TT", "TTxGT>AC", "TTxGT>AG", "TTxGT>AT", "TTxGT>CG", "TTxGT>CT", "TTxGT>GT"},
    {"ACxAA>AA", "ACxAA>CC", "ACxAA>GG", "ACxAA>TT", "ACxAA>AC", "ACxAA>AG", "ACxAA>AT", "ACxAA>CG", "ACxAA>CT", "ACxAA>GT"},
    {"ACxCC>AA", "ACxCC>CC", "ACxCC>GG", "ACxCC>TT", "ACxCC>AC", "ACxCC>AG", "ACxCC>AT", "ACxCC>CG", "ACxCC>CT", "ACxCC>GT"},
    {"ACxGG>AA", "ACxGG>CC", "ACxGG>GG", "ACxGG>TT", "ACxGG>AC", "ACxGG>AG", "ACxGG>AT", "ACxGG>CG", "ACxGG>CT", "ACxGG>GT"},
    {"ACxTT>AA", "ACxTT>CC", "ACxTT>GG", "ACxTT>TT", "ACxTT>AC", "ACxTT>AG", "ACxTT>AT", "ACxTT>CG", "ACxTT>CT", "ACxTT>GT"},
    {"ACxAC>AA", "ACxAC>CC", "ACxAC>GG", "ACxAC>TT", "ACxAC>AC", "ACxAC>AG", "ACxAC>AT", "ACxAC>CG", "ACxAC>CT", "ACxAC>GT"},
    {"ACxAG>AA", "ACxAG>CC", "ACxAG>GG", "ACxAG>TT", "ACxAG>AC", "ACxAG>AG", "ACxAG>AT", "ACxAG>CG", "ACxAG>CT", "ACxAG>GT"},
    {"ACxAT>AA", "ACxAT>CC", "ACxAT>GG", "ACxAT>TT", "ACxAT>AC", "ACxAT>AG", "ACxAT>AT", "ACxAT>CG", "ACxAT>CT", "ACxAT>GT"},
    {"ACxCG>AA", "ACxCG>CC", "ACxCG>GG", "ACxCG>TT", "ACxCG>AC", "ACxCG>AG", "ACxCG>AT", "ACxCG>CG", "ACxCG>CT", "ACxCG>GT"},
    {"ACxCT>AA", "ACxCT>CC", "ACxCT>GG", "ACxCT>TT", "ACxCT>AC", "ACxCT>AG", "ACxCT>AT", "ACxCT>CG", "ACxCT>CT", "ACxCT>GT"},
    {"ACxGT>AA", "ACxGT>CC", "ACxGT>GG", "ACxGT>TT", "ACxGT>AC", "ACxGT>AG", "ACxGT>AT", "ACxGT>CG", "ACxGT>CT", "ACxGT>GT"},
    {"AGxAA>AA", "AGxAA>CC", "AGxAA>GG", "AGxAA>TT", "AGxAA>AC", "AGxAA>AG", "AGxAA>AT", "AGxAA>CG", "AGxAA>CT", "AGxAA>GT"},
    {"AGxCC>AA", "AGxCC>CC", "AGxCC>GG", "AGxCC>TT", "AGxCC>AC", "AGxCC>AG", "AGxCC>AT", "AGxCC>CG", "AGxCC>CT", "AGxCC>GT"},
    {"AGxGG>AA", "AGxGG>CC", "AGxGG>GG", "AGxGG>TT", "AGxGG>AC", "AGxGG>AG", "AGxGG>AT", "AGxGG>CG", "AGxGG>CT", "AGxGG>GT"},
    {"AGxTT>AA", "AGxTT>CC", "AGxTT>GG", "AGxTT>TT", "AGxTT>AC", "AGxTT>AG", "AGxTT>AT", "AGxTT>CG", "AGxTT>CT", "AGxTT>GT"},
    {"AGxAC>AA", "AGxAC>CC", "AGxAC>GG", "AGxAC>TT", "AGxAC>AC", "AGxAC>AG", "AGxAC>AT", "AGxAC>CG", "AGxAC>CT", "AGxAC>GT"},
    {"AGxAG>AA", "AGxAG>CC", "AGxAG>GG", "AGxAG>TT", "AGxAG>AC", "AGxAG>AG", "AGxAG>AT", "AGxAG>CG", "AGxAG>CT", "AGxAG>GT"},
    {"AGxAT>AA", "AGxAT>CC", "AGxAT>GG", "AGxAT>TT", "AGxAT>AC", "AGxAT>AG", "AGxAT>AT", "AGxAT>CG", "AGxAT>CT", "AGxAT>GT"},
    {"AGxCG>AA", "AGxCG>CC", "AGxCG>GG", "AGxCG>TT", "AGxCG>AC", "AGxCG>AG", "AGxCG>AT", "AGxCG>CG", "AGxCG>CT", "AGxCG>GT"},
    {"AGxCT>AA", "AGxCT>CC", "AGxCT>GG", "AGxCT>TT", "AGxCT>AC", "AGxCT>AG", "AGxCT>AT", "AGxCT>CG", "AGxCT>CT", "AGxCT>GT"},
    {"AGxGT>AA", "AGxGT>CC", "AGxGT>GG", "AGxGT>TT", "AGxGT>AC", "AGxGT>AG", "AGxGT>AT", "AGxGT>CG", "AGxGT>CT", "AGxGT>GT"},
    {"ATxAA>AA", "ATxAA>CC", "ATxAA>GG", "ATxAA>TT", "ATxAA>AC", "ATxAA>AG", "ATxAA>AT", "ATxAA>CG", "ATxAA>CT", "ATxAA>GT"},
    {"ATxCC>AA", "ATxCC>CC", "ATxCC>GG", "ATxCC>TT", "ATxCC>AC", "ATxCC>AG", "ATxCC>AT", "ATxCC>CG", "ATxCC>CT", "ATxCC>GT"},
    {"ATxGG>AA", "ATxGG>CC", "ATxGG>GG", "ATxGG>TT", "ATxGG>AC", "ATxGG>AG", "ATxGG>AT", "ATxGG>CG", "ATxGG>CT", "ATxGG>GT"},
    {"ATxTT>AA", "ATxTT>CC", "ATxTT>GG", "ATxTT>TT", "ATxTT>AC", "ATxTT>AG", "ATxTT>AT", "ATxTT>CG", "ATxTT>CT", "ATxTT>GT"},
    {"ATxAC>AA", "ATxAC>CC", "ATxAC>GG", "ATxAC>TT", "ATxAC>AC", "ATxAC>AG", "ATxAC>AT", "ATxAC>CG", "ATxAC>CT", "ATxAC>GT"},
    {"ATxAG>AA", "ATxAG>CC", "ATxAG>GG", "ATxAG>TT", "ATxAG>AC", "ATxAG>AG", "ATxAG>AT", "ATxAG>CG", "ATxAG>CT", "ATxAG>GT"},
    {"ATxAT>AA", "ATxAT>CC", "ATxAT>GG", "ATxAT>TT", "ATxAT>AC", "ATxAT>AG", "ATxAT>AT", "ATxAT>CG", "ATxAT>CT", "ATxAT>GT"},
    {"ATxCG>AA", "ATxCG>CC", "ATxCG>GG", "ATxCG>TT", "ATxCG>AC", "ATxCG>AG", "ATxCG>AT", "ATxCG>CG", "ATxCG>CT", "ATxCG>GT"},
    {"ATxCT>AA", "ATxCT>CC", "ATxCT>GG", "ATxCT>TT", "ATxCT>AC", "ATxCT>AG", "ATxCT>AT", "ATxCT>CG", "ATxCT>CT", "ATxCT>GT"},
    {"ATxGT>AA", "ATxGT>CC", "ATxGT>GG", "ATxGT>TT", "ATxGT>AC", "ATxGT>AG", "ATxGT>AT", "ATxGT>CG", "ATxGT>CT", "ATxGT>GT"},
    {"CGxAA>AA", "CGxAA>CC", "CGxAA>GG", "CGxAA>TT", "CGxAA>AC", "CGxAA>AG", "CGxAA>AT", "CGxAA>CG", "CGxAA>CT", "CGxAA>GT"},
    {"CGxCC>AA", "CGxCC>CC", "CGxCC>GG", "CGxCC>TT", "CGxCC>AC", "CGxCC>AG", "CGxCC>AT", "CGxCC>CG", "CGxCC>CT", "CGxCC>GT"},
    {"CGxGG>AA", "CGxGG>CC", "CGxGG>GG", "CGxGG>TT", "CGxGG>AC", "CGxGG>AG", "CGxGG>AT", "CGxGG>CG", "CGxGG>CT", "CGxGG>GT"},
    {"CGxTT>AA", "CGxTT>CC", "CGxTT>GG", "CGxTT>TT", "CGxTT>AC", "CGxTT>AG", "CGxTT>AT", "CGxTT>CG", "CGxTT>CT", "CGxTT>GT"},
    {"CGxAC>AA", "CGxAC>CC", "CGxAC>GG", "CGxAC>TT", "CGxAC>AC", "CGxAC>AG", "CGxAC>AT", "CGxAC>CG", "CGxAC>CT", "CGxAC>GT"},
    {"CGxAG>AA", "CGxAG>CC", "CGxAG>GG", "CGxAG>TT", "CGxAG>AC", "CGxAG>AG", "CGxAG>AT", "CGxAG>CG", "CGxAG>CT", "CGxAG>GT"},
    {"CGxAT>AA", "CGxAT>CC", "CGxAT>GG", "CGxAT>TT", "CGxAT>AC", "CGxAT>AG", "CGxAT>AT", "CGxAT>CG", "CGxAT>CT", "CGxAT>GT"},
    {"CGxCG>AA", "CGxCG>CC", "CGxCG>GG", "CGxCG>TT", "CGxCG>AC", "CGxCG>AG", "CGxCG>AT", "CGxCG>CG", "CGxCG>CT", "CGxCG>GT"},
    {"CGxCT>AA", "CGxCT>CC", "CGxCT>GG", "CGxCT>TT", "CGxCT>AC", "CGxCT>AG", "CGxCT>AT", "CGxCT>CG", "CGxCT>CT", "CGxCT>GT"},
    {"CGxGT>AA", "CGxGT>CC", "CGxGT>GG", "CGxGT>TT", "CGxGT>AC", "CGxGT>AG", "CGxGT>AT", "CGxGT>CG", "CGxGT>CT", "CGxGT>GT"},
    {"CTxAA>AA", "CTxAA>CC", "CTxAA>GG", "CTxAA>TT", "CTxAA>AC", "CTxAA>AG", "CTxAA>AT", "CTxAA>CG", "CTxAA>CT", "CTxAA>GT"},
    {"CTxCC>AA", "CTxCC>CC", "CTxCC>GG", "CTxCC>TT", "CTxCC>AC", "CTxCC>AG", "CTxCC>AT", "CTxCC>CG", "CTxCC>CT", "CTxCC>GT"},
    {"CTxGG>AA", "CTxGG>CC", "CTxGG>GG", "CTxGG>TT", "CTxGG>AC", "CTxGG>AG", "CTxGG>AT", "CTxGG>CG", "CTxGG>CT", "CTxGG>GT"},
    {"CTxTT>AA", "CTxTT>CC", "CTxTT>GG", "CTxTT>TT", "CTxTT>AC", "CTxTT>AG", "CTxTT>AT", "CTxTT>CG", "CTxTT>CT", "CTxTT>GT"},
    {"CTxAC>AA", "CTxAC>CC", "CTxAC>GG", "CTxAC>TT", "CTxAC>AC", "CTxAC>AG", "CTxAC>AT", "CTxAC>CG", "CTxAC>CT", "CTxAC>GT"},
    {"CTxAG>AA", "CTxAG>CC", "CTxAG>GG", "CTxAG>TT", "CTxAG>AC", "CTxAG>AG", "CTxAG>AT", "CTxAG>CG", "CTxAG>CT", "CTxAG>GT"},
    {"CTxAT>AA", "CTxAT>CC", "CTxAT>GG", "CTxAT>TT", "CTxAT>AC", "CTxAT>AG", "CTxAT>AT", "CTxAT>CG", "CTxAT>CT", "CTxAT>GT"},
    {"CTxCG>AA", "CTxCG>CC", "CTxCG>GG", "CTxCG>TT", "CTxCG>AC", "CTxCG>AG", "CTxCG>AT", "CTxCG>CG", "CTxCG>CT", "CTxCG>GT"},
    {"CTxCT>AA", "CTxCT>CC", "CTxCT>GG", "CTxCT>TT", "CTxCT>AC", "CTxCT>AG", "CTxCT>AT", "CTxCT>CG", "CTxCT>CT", "CTxCT>GT"},
    {"CTxGT>AA", "CTxGT>CC", "CTxGT>GG", "CTxGT>TT", "CTxGT>AC", "CTxGT>AG", "CTxGT>AT", "CTxGT>CG", "CTxGT>CT", "CTxGT>GT"},
    {"GTxAA>AA", "GTxAA>CC", "GTxAA>GG", "GTxAA>TT", "GTxAA>AC", "GTxAA>AG", "GTxAA>AT", "GTxAA>CG", "GTxAA>CT", "GTxAA>GT"},
    {"GTxCC>AA", "GTxCC>CC", "GTxCC>GG", "GTxCC>TT", "GTxCC>AC", "GTxCC>AG", "GTxCC>AT", "GTxCC>CG", "GTxCC>CT", "GTxCC>GT"},
    {"GTxGG>AA", "GTxGG>CC", "GTxGG>GG", "GTxGG>TT", "GTxGG>AC", "GTxGG>AG", "GTxGG>AT", "GTxGG>CG", "GTxGG>CT", "GTxGG>GT"},
    {"GTxTT>AA", "GTxTT>CC", "GTxTT>GG", "GTxTT>TT", "GTxTT>AC", "GTxTT>AG", "GTxTT>AT", "GTxTT>CG", "GTxTT>CT", "GTxTT>GT"},
    {"GTxAC>AA", "GTxAC>CC", "GTxAC>GG", "GTxAC>TT", "GTxAC>AC", "GTxAC>AG", "GTxAC>AT", "GTxAC>CG", "GTxAC>CT", "GTxAC>GT"},
    {"GTxAG>AA", "GTxAG>CC", "GTxAG>GG", "GTxAG>TT", "GTxAG>AC", "GTxAG>AG", "GTxAG>AT", "GTxAG>CG", "GTxAG>CT", "GTxAG>GT"},
    {"GTxAT>AA", "GTxAT>CC", "GTxAT>GG", "GTxAT>TT", "GTxAT>AC", "GTxAT>AG", "GTxAT>AT", "GTxAT>CG", "GTxAT>CT", "GTxAT>GT"},
    {"GTxCG>AA", "GTxCG>CC", "GTxCG>GG", "GTxCG>TT", "GTxCG>AC", "GTxCG>AG", "GTxCG>AT", "GTxCG>CG", "GTxCG>CT", "GTxCG>GT"},
    {"GTxCT>AA", "GTxCT>CC", "GTxCT>GG", "GTxCT>TT", "GTxCT>AC", "GTxCT>AG", "GTxCT>AT", "GTxCT>CG", "GTxCT>CT", "GTxCT>GT"},
    {"GTxGT>AA", "GTxGT>CC", "GTxGT>GG", "GTxGT>TT", "GTxGT>AC", "GTxGT>AG", "GTxGT>AT", "GTxGT>CG", "GTxGT>CT", "GTxGT>GT"}
};

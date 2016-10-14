/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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


#define BOOST_TEST_MODULE dng::relationship_graph_peeling

#include <dng/relationship_graph.h>

#include "../boost_test_helper.h"
#include "fixture_trio_workspace.h"
#include "relationship_graph_helper.h"


//struct FixturePedigree : public ReadTrioFromFile{
//
//
//    dng::RelationshipGraph relationship_graph;
//
//    FixturePedigree(std::string s = "FixturePedigree") : ReadTrioFromFile(s) {
//        BOOST_TEST_MESSAGE("set up fixture: " << fixture);
//        relationship_graph.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
//    }
//
//    ~FixturePedigree() {
//        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
//    }
//};


namespace dng {
BOOST_FIXTURE_TEST_CASE(test_peeling_forward_each_op, TrioWorkspace ) {

    TransitionVector mat;
    mat.assign(5, {});
    auto f81_matrix_3 = f81::matrix(3e-8, {{0.3,0.2,0.2,0.3}});
    auto f81_matrix_2 = f81::matrix(2e-8, {{0.3,0.2,0.2,0.3}});

    GenotypeArray expected_lower_mom(10);
    GenotypeArray expected_lower_dad(10);
    GenotypeArray expected_lower_root(10);
    GenotypeArray expected_root_product(10);
    double expected_ret = 0;


    int num_mutation = -1;
    expected_lower_mom << 8.778182862e-16, 8.778182862e-16, 8.388581887e-08, 1.168700868e-13, 8.778182862e-16, 8.388581887e-08, 1.168700868e-13, 0.9999999568, 2.153347088e-05, 2.328623553e-13;
    expected_lower_dad << 1.2241233419508e-08, 2.50389467694902e-08, 1.10660049747288e-07, 0.499400615090229, 3.78366601194724e-08, 1.2345776309727e-07, 0.499400627887942, 2.09078866075068e-07, 0.499400713509045, 0.998801217939224;
    expected_lower_root << 3.05392180906051e-23, 6.24667327177925e-23, 2.58056140979112e-14, 1.24589604463533e-15, 9.43942473449798e-23, 2.87900050574266e-14, 1.24589607656285e-15, 2.09078857033822e-07, 1.16458849625194e-07, 2.49179205873145e-15;
    expected_root_product << 4.57539160330409e-27, 3.74238850910963e-30, 7.73163125352332e-18, 1.11962644069753e-22, 9.42716715512861e-27, 5.75052395049989e-18, 7.46417646259496e-23, 2.08828238101239e-07, 3.48922865425786e-11, 3.73320771634523e-19;
    expected_ret = -15.38158668;

    mat[2] = meiosis_diploid_matrix(f81_matrix_3, f81_matrix_3, num_mutation);
    mat[3] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    mat[4] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    workspace.CleanupFast();
    (*r_graph.peeling_functions_[0])(workspace, r_graph.family_members_[0], mat);
    boost_check_matrix(expected_lower_mom, workspace.lower[1]);
    (*r_graph.peeling_functions_[1])(workspace, r_graph.family_members_[1], mat);
    boost_check_matrix(expected_lower_dad, workspace.lower[0]);
    (*r_graph.peeling_functions_[2])(workspace, r_graph.family_members_[2],mat);
    boost_check_matrix(expected_lower_root, workspace.lower[0]);
    GenotypeArray root_product = workspace.lower[0] * workspace.upper[0];
    boost_check_matrix(expected_root_product, root_product);
    auto ret = log(root_product.sum());
    BOOST_CHECK_CLOSE(expected_ret, ret, BOOST_CLOSE_PERCENTAGE_THRESHOLD);


    num_mutation = 0;
    expected_lower_mom << 1.62968750071628e-19, 1.62968749631172e-19, 7.84802386914498e-08, 4.87276562714167e-17, 1.62968749190716e-19, 7.84802384793411e-08, 4.87276561397203e-17, 0.999999956756758, 2.15280654759117e-05, 9.7292343792761e-17;
    expected_lower_dad << 9.36491392262008e-11, 1.28913624242886e-08, 9.24382700267887e-08, 0.499400609015176, 2.5689075709351e-08, 1.05235983311851e-07, 0.499400621812889, 1.84782890914351e-07, 0.499400701359796, 0.998801217936702;
    expected_lower_root << 2.76293563891165e-28, 3.80334564462965e-26, 2.10566794406379e-14, 1.47338432808165e-18, 7.5790619124609e-26, 2.39718934548084e-14, 1.4733843618567e-18, 1.8478288292374e-07, 1.13759382103489e-07, 2.94676865588702e-18;
    expected_root_product << 4.13943555635268e-32, 2.2785883649362e-33, 6.30880087724149e-18, 1.32405914452713e-25, 7.56922011006824e-30, 4.78815294320859e-18, 8.82706116586081e-26, 1.84561387122897e-07, 3.40834979053659e-11, 4.41485454048767e-22;
    expected_ret = -15.50509905;

    mat[2] = meiosis_diploid_matrix(f81_matrix_3, f81_matrix_3, num_mutation);
    mat[3] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    mat[4] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    workspace.CleanupFast();
    (*r_graph.peeling_functions_[0])(workspace, r_graph.family_members_[0], mat);
    boost_check_matrix(expected_lower_mom, workspace.lower[1]);
    (*r_graph.peeling_functions_[1])(workspace, r_graph.family_members_[1], mat);
    boost_check_matrix(expected_lower_dad, workspace.lower[0]);
    (*r_graph.peeling_functions_[2])(workspace, r_graph.family_members_[2],mat);
    boost_check_matrix(expected_lower_root, workspace.lower[0]);
    root_product = workspace.lower[0] * workspace.upper[0];
    boost_check_matrix(expected_root_product, root_product);
    ret = log(root_product.sum());
    BOOST_CHECK_CLOSE(expected_ret, ret, BOOST_CLOSE_PERCENTAGE_THRESHOLD);



    num_mutation = 1;
    expected_lower_mom << 8.48435019136673e-16, 8.4843501799058e-16, 5.4055801825711e-09, 1.16792139751089e-13, 8.48435016844487e-16, 5.40558039467985e-09, 1.16792139750074e-13, 3.51224853319191e-13, 5.40540626701447e-09, 2.3273584448304e-13;
    expected_lower_dad << 3.70953848566537e-20, 3.70953853807702e-20, 5.72844550473873e-20, 3.05851180397549e-20, 3.70953859048867e-20, 5.72844555715038e-20, 3.05851185638714e-20, 7.74735252381209e-20, 5.07741882304885e-20, 2.40748512228561e-20;
    expected_lower_root << 9.1351614607608e-35, 9.13516157751524e-35, 3.09645867899861e-28, 7.53193403085758e-35, 9.13516169426968e-35, 3.096459060003e-28, 7.53193414977335e-35, 7.63149553815086e-34, 2.7445521771194e-28, 5.92870660095436e-35;
    expected_root_product << 1.36863166955968e-38, 5.47288488274259e-42, 9.27731330359118e-32, 6.76858436692924e-42, 9.12329921613473e-39, 6.18487629674674e-32, 4.51238964919539e-42, 7.62234780655847e-34, 8.22296470412592e-32, 8.88239977853428e-39;
    expected_ret = -70.51466136;

    mat[2] = meiosis_diploid_matrix(f81_matrix_3, f81_matrix_3, num_mutation);
    mat[3] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    mat[4] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    workspace.CleanupFast();
    (*r_graph.peeling_functions_[0])(workspace, r_graph.family_members_[0], mat);
    boost_check_matrix(expected_lower_mom, workspace.lower[1]);
    (*r_graph.peeling_functions_[1])(workspace, r_graph.family_members_[1], mat);
    boost_check_matrix(expected_lower_dad, workspace.lower[0]);
    (*r_graph.peeling_functions_[2])(workspace, r_graph.family_members_[2],mat);
    boost_check_matrix(expected_lower_root, workspace.lower[0]);
    root_product = workspace.lower[0] * workspace.upper[0];
    boost_check_matrix(expected_root_product, root_product);
    ret = log(root_product.sum());
    BOOST_CHECK_CLOSE(expected_ret, ret, BOOST_CLOSE_PERCENTAGE_THRESHOLD);



    num_mutation = 2;
    expected_lower_mom << 2.92202983431692e-17, 2.9220299489703e-17, 9.49256366449536e-22, 2.92193582590732e-17, 2.92203006362368e-17, 9.49256366450609e-22, 2.92193594056071e-17, 1.71039811434089e-32, 9.49256366441022e-22, 2.92184181749773e-17;
    expected_lower_dad << 1.87172063020916e-36, 1.87172061042194e-36, 1.51175764174171e-36, 1.29582332367199e-36, 1.87172059063473e-36, 1.51175762195449e-36, 1.29582330388478e-36, 1.15179465327426e-36, 9.35860335204542e-37, 7.19926017134825e-37;
    expected_lower_root << 5.46887569084554e-53, 5.46887625591219e-53, 4.02473034425802e-59, 3.78619360180291e-53, 5.46887682097884e-53, 4.02473029160487e-59, 3.78619397521919e-53, 1.58861465907825e-69, 2.49152733552297e-59, 2.10351151276028e-53;
    expected_root_product << 8.19348021326919e-57, 3.27640950109093e-60, 1.20585120736142e-62, 3.4024688371325e-60, 5.46177520265145e-57, 8.03900794390654e-63, 2.26831278180242e-60, 1.58671041626866e-69, 7.4648758767172e-63, 3.15148504600963e-57;
    expected_ret = -128.4250364;

    mat[2] = meiosis_diploid_matrix(f81_matrix_3, f81_matrix_3, num_mutation);
    mat[3] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    mat[4] = mitosis_diploid_matrix(f81_matrix_2, num_mutation);
    workspace.CleanupFast();
    (*r_graph.peeling_functions_[0])(workspace, r_graph.family_members_[0], mat);
    boost_check_matrix(expected_lower_mom, workspace.lower[1]);
    (*r_graph.peeling_functions_[1])(workspace, r_graph.family_members_[1], mat);
    boost_check_matrix(expected_lower_dad, workspace.lower[0]);
    (*r_graph.peeling_functions_[2])(workspace, r_graph.family_members_[2],mat);
    boost_check_matrix(expected_lower_root, workspace.lower[0]);
    root_product = workspace.lower[0] * workspace.upper[0];
    boost_check_matrix(expected_root_product, root_product);
    ret = log(root_product.sum());
    BOOST_CHECK_CLOSE(expected_ret, ret, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

}
} // namespace dng


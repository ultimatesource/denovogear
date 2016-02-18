//
// Created by steven on 1/24/16.
//

#ifndef DENOVOGEAR_TESTH_HELPER_PEELING_H
#define DENOVOGEAR_TESTH_HELPER_PEELING_H

#include <dng/peeling.h>

//TODO: Fix this files!! New filename system and mixed global and local here.
std::random_device rd;
std::mt19937 random_gen_mt(rd());

struct Fx {

    const int CHILD_OFFSET = 2;
    std::string fixture;

    dng::peel::family_members_t family;
    std::vector<dng::TransitionMatrix> trans_matrix;
    std::vector<dng::GenotypeArray> upper_array;
    std::vector<dng::GenotypeArray> lower_array;

    int num_child;
    int total_family_size;
    std::uniform_int_distribution<> rand_unif;

    Fx(std::string s = "") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture " << s);
        rand_unif = std::uniform_int_distribution<>(1,10);
        init_family();
    }

    void init_family(){

        int num_child = rand_unif(random_gen_mt);
        total_family_size = num_child + CHILD_OFFSET;

        family.clear();
        trans_matrix.resize(total_family_size);
        upper_array.resize(total_family_size);
        lower_array.resize(total_family_size);

        for (int k = 0; k < total_family_size; ++k) {
            family.push_back(k);
            lower_array[k] = dng::GenotypeArray::Random();
            if (k < CHILD_OFFSET) {
                trans_matrix[k] = dng::TransitionMatrix::Random(10, 10);
                upper_array[k] = dng::GenotypeArray::Random();
            }
            else {
                trans_matrix[k] = dng::TransitionMatrix::Random(100, 10);
            }

        }
    }

    void init_family_parent_child_only(){

        total_family_size = 2;

        family.clear();
        trans_matrix.resize(total_family_size);
        upper_array.resize(total_family_size);
        lower_array.resize(total_family_size);

        for (int k = 0; k < total_family_size; ++k) {
            family.push_back(k);
            trans_matrix[k] = dng::TransitionMatrix::Random(10, 10);
            lower_array[k] = dng::GenotypeArray::Random();
//            upper_array[k] = GenotypeArray::Random();

        }
    }

    ~Fx() {
        BOOST_TEST_MESSAGE("tear down fixture " << fixture);
    }
};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }




void copy_family_to_workspace(dng::peel::workspace_t &workspace, dng::TransitionVector &full_matrix,
                              int total_family_size, const std::vector<dng::GenotypeArray> &lower,
                              const std::vector<dng::GenotypeArray> &upper,
                              const std::vector<dng::TransitionMatrix> &trans_mat) {
    workspace.Resize(total_family_size);
    full_matrix.resize(total_family_size);
    for (int l = 0; l < total_family_size; ++l) {
        workspace.lower[l] = lower[l];
        workspace.upper[l] = upper[l];
        full_matrix[l] = trans_mat[l];
    }
}

#endif //DENOVOGEAR_TESTH_HELPER_PEELING_H

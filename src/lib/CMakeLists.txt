add_library(libdng STATIC
  ad.cc
  find_mutations.cc
  inheritance_model.cc
  likelihood.cc
  mutation_stats.cc
  mutation.cc
  newick.cc
  peeling.cc
  regions.cc
  relationship_graph.cc
  read_group.cc
  stats.cc
  utility.cc
  task/call.cc
  task/loglike.cc
  task/pileup.cc
)

target_include_directories(libdng PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_include_directories(libdng PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/..")
target_link_libraries(libdng
  HTSLIB::HTSLIB
  Threads::Threads
  EIGEN3::EIGEN3
  Boost::PROGRAM_OPTIONS
  Boost::FILESYSTEM
  Boost::SYSTEM
)

set_target_properties(libdng PROPERTIES OUTPUT_NAME dng)

if(TARGET ext_eigen3)
  add_dependencies(libdng ext_eigen3)
endif()

# Some compilers have incomplete c++11 support provide a list of
# common features to toss an error if someone tries to use
# an early compiler
target_compile_features(libdng PUBLIC
  cxx_auto_type
  cxx_constexpr
  cxx_decltype
  cxx_deleted_functions
  cxx_lambdas
  cxx_nonstatic_member_init
  cxx_nullptr
  cxx_range_for
  cxx_rvalue_references
  cxx_strong_enums
  cxx_trailing_return_types
  cxx_attributes
  cxx_reference_qualified_functions
  cxx_decltype_incomplete_return_types
)

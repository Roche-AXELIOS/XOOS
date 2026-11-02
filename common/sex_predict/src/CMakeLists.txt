add_default_library(xoos_sex_predict
        xoos/sex_predict/allele-counter.cpp
        xoos/sex_predict/cascade-contig-span.cpp
        xoos/sex_predict/cascade-depth-profile.cpp
        xoos/sex_predict/cascade-heterozygosity.cpp
        xoos/sex_predict/cascade-residual-chry.cpp
        xoos/sex_predict/cascade-resolve.cpp
        xoos/sex_predict/scan-chunk.cpp
        xoos/sex_predict/sex.cpp
        xoos/sex_predict/sex-predictor.cpp
        xoos/sex_predict/sex-predictor-cascade.cpp
        xoos/sex_predict/vcf-sex-predictor.cpp
)

target_include_directories(xoos_sex_predict PRIVATE .)
target_include_directories(xoos_sex_predict SYSTEM INTERFACE .)

find_package(xoos-htslib REQUIRED)
find_package(xoos-taskflow REQUIRED)
target_link_libraries(xoos_sex_predict
        PUBLIC
        xoos-htslib::xoos-htslib
        PRIVATE
        xoos-taskflow::xoos-taskflow
)

target_link_xoos_libraries(xoos_sex_predict
        PUBLIC
        io
        PRIVATE
        error
        genome
        log
        util
)

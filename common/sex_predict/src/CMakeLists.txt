add_library(xoos_sex_predict
        STATIC
        xoos/sex_predict/sex.cpp
        xoos/sex_predict/sex-predictor.cpp
        xoos/sex_predict/vcf-sex-predictor.cpp
)

target_include_directories(xoos_sex_predict PRIVATE .)
target_include_directories(xoos_sex_predict SYSTEM INTERFACE .)

find_package(xoos-htslib REQUIRED)
target_link_libraries(xoos_sex_predict
        PUBLIC
        xoos-htslib::xoos-htslib
        PRIVATE
)

target_link_xoos_libraries(xoos_sex_predict
        PUBLIC
        io
        PRIVATE
        error
        log
        util
)

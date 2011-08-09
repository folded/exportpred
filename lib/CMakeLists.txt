project(ghmm)

set(ghmm_VERSION_MAJOR 2)
set(ghmm_VERSION_MINOR 0)
set(ghmm_VERSION_PATCH 0)

set(GHMM_VERSION ${ghmm_VERSION_MAJOR}.${ghmm_VERSION_MINOR}.${ghmm_VERSION_PATCH})

option(GHMM_DEBUG "Compile libghmm in debug code" OFF)

include_directories("${exportpred_SOURCE_DIR}/include")

add_library(ghmm
            ghmm.cc
            ghmm_util.cc
           )

set_target_properties(ghmm PROPERTIES
                      VERSION   "${ghmm_VERSION_MAJOR}.${ghmm_VERSION_MINOR}.${ghmm_VERSION_PATCH}"
                      SOVERSION "${ghmm_VERSION_MAJOR}.${ghmm_VERSION_MINOR}")

install(TARGETS ghmm
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

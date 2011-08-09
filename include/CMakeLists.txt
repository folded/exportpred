install(DIRECTORY GHMM
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
        FILES_MATCHING
          PATTERN "*.hh"
          PATTERN "internal" EXCLUDE
)

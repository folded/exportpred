include_directories("${exportpred_SOURCE_DIR}/include")

add_executable       (exportpred exportpred.cc ss_model.cc signalp_model.cc)
target_link_libraries(exportpred ghmm)

foreach(tgt exportpred)
  install(TARGETS ${tgt}
          RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
endforeach(tgt)

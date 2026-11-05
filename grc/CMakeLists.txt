install(FILES bonded_usrp_bonded_source.block.yml
              bonded_usrp_overlap_reconstructor.block.yml
              bonded_usrp_bonded_sink.block.yml
              bonded_usrp_band_splitter.block.yml
              bonded_usrp_ofdm_splitter.block.yml
    DESTINATION ${GR_PKG_DATA_DIR}/grc/blocks)

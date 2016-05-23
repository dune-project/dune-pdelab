# TODO add documentation
function(add_dune_vcl_flags _targets)
  if(VCL_FOUND)
    foreach(_target ${_targets})
      target_include_directories(${_target} PUBLIC ${VCL_INCLUDE_DIR})
    endforeach(_target ${_targets})
  end(VCL_FOUND)
endfunction(add_dune_vcl_flags)

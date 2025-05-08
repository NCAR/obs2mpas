set( obs2mpas_src
	setvar.F90 
	netcdf_mod.f90
	goes_abi_mod.F90
	himawari_ahi_mod.F90
	read_write_indx_mod.F90
	read_write_mpas_mod.F90
	utils_mod.F90
	main.F90
)

add_executable(obs2mpas.x  ${obs2mpas_src} )

set(FORTRAN_COMPILER_GNU_FLAGS
    $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none -mcmodel=medium>)
target_compile_options(obs2mpas.x PUBLIC ${FORTRAN_COMPILER_GNU_FLAGS})

target_link_libraries(obs2mpas.x
  PUBLIC
  NetCDF::NetCDF_Fortran
  atlas_f
)

## v2.4.0
- Integrated with Sectral Transform v2 (mkv, anandogc)

### Bug fixes
- Resolved(#2): In para.yaml, io.time.*save_first are now relative to time.init (siddhesh, anandogc)
- Resolved(#1): HDF5 IO now gives error and the program is terminated on any IO error (biplab_dutta, anandogc)
- output_ET.cc: shell_to_shell_file flushed properly. (ksreddy)


## v2.3.1

### Added features:
- STRATIFIED added as a program kind. (abhishek_kir)
- Necessary functions for GPU are declared, definitions will be added later. (mkv, anandogc)

### Bug fixes
- compute_force_random.cc - FORCE::Compute_force_using_random_energy_helicity_spectrum, FORCE::Compute_force_using_random_energy_helicity_spectrum line no 241 and 305. (abhishek_kir)
- void FORCE::Compute_force_RBC_basic_add(FluidVF& U, FluidSF& T) line 166: global.program.basis_type == "ChFF" section completely changed. (mkv)
- void FORCE::Compute_force_RBC_basic_assign(FluidVF& U, FluidSF& T): /TWO (mkv)
- DP Correlation::Get_Nusselt_no(FluidVF& U, FluidSF& T) factor 2 to 4 in line 163 and 168 (mkv)
- Global.cc: energy_transfer.flux.no_spheres + 1 (anandogc)
- void FORCE::Compute_force_stratified_random(FluidVF& U, FluidSF& T) (abhishek_kir)
- ISclar_main.cc DIAGNOSTIC Real outout (biplab)
- ISclar_main.cc (global.program.basis_type == "SFF") (mkv)
- GROSSMANN_LOHSE: dt bug fixed (anandogc, pambr)
- GROSSMANN_LOHSE: Only master was writing. Now a MPI_Reduce is done before writing. (anandogc, pambr)

### Revisions:
   - Earlier even if init condition was not from mode it parsed modes, now it doesn't. (anandogc)
   - Earlier even if init condition was not from rediced it parsed N_in_reduced and N_out_reduced, now it doesn't. (anandogc)
   - GROSSMANN_LOHSE: Comment in output file (misc.d) changed. (anandogc, pambr)
   - HDF5 Reader and Writer generalized so that 'dataset name' can now be different from 'file name' (anandogc)

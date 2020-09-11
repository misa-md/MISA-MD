<a name="unreleased"></a>
## [Unreleased]


<a name="v0.4.0"></a>
## [v0.4.0] - 2020-09-11
### Ci
- **gitlab-ci:** setup gitlab-ci for pkg version 0.4.x

### Docs
- **changelog:** update changelog for v0.4.0

### Feat
- **arch:** add a new arch func ${ARCH_NAME}_nei_offset_init to init neighbor offset on device side

### Refactor
- **arch:** change func signature of archAccDomainInit, and move its call to class simulation
- **arch:** rename signature of architecture acceleration functions that called at host side
- **arch:** pass higher-level potential and domain parameters into arch acceleration api
- **neighbor:** make neighbor offset range radius configurable

### BREAKING CHANGE

the arch of initializing domain has changed: ${ARCH_NAME}_accelerate_init -> ${ARCH_NAME}_domain_init.

the arch api has changed.


<a name="v0.4.0-beta.4"></a>
## [v0.4.0-beta.4] - 2020-05-21
### Build
- **cmake:** add `MD_` prefix to variable names of cmake options
- **pkg:** upgrade pkg.yaml file to version 2

### Docs
- **changelog:** changlog for hotfix of "wrong parameters passing to func simulation::createAtoms"
- **changelog:** add changelog for hotfix of `MD_` prefix of cmake variables
- **changelog:** update changelog for v0.4.0
- **changelog:** add commit types of "build" and "improvement" to CHANGELOG.md file

### Fix
- **cmake:** add missing `_MD` prefix to cmake variabls in arch/arch_configure.cmake file
- **config:** fix empty value of term `creation.create_t_set` in example/config.yaml file
- **simulation:** fix wrong parameters passing (tset and create seed) to func simulation::createAtoms

### Improvement
- **eam:** remove unnecessary temporary variables in eam computing
- **logs:** log execution time of the whole program, not only time of timestep loop
- **neighbor:** ues 1d linear offset(instead of 3d offset) in neighbor inde and its iterator

### Refactor
- **eam:** remove unnecessary BccLattice::IndexOf3DIndex calling while computing eam
- **neighbor:** simplify neighbor index calculation and compare it with old method in unit test
- **neighbor:** make members in NeiIterator constant and adjust comparison order in == operator

### BREAKING CHANGE

when performing cmake configuration, the option names have `MD_` prefix now
(e.g. change 'SUNWAY_ARCH_ENABLE_FLAG' to `MD_SUNWAY_ARCH_ENABLE_FLAG`).


<a name="v0.4.0-beta.3"></a>
## [v0.4.0-beta.3] - 2019-11-26
### Build
- **arch:** move architecture related files and cmake config to src/arch directory
- **pkg:** update libcomm to version 0.3.2
- **pkg:** update libcomm to version 0.3.1
- **pkg:** update googletest to version 1.10.0

### Ci
- **gitlab-ci:** update gitlab ci cache key due to dependency changed

### Docs
- **changelog:** update changelog for v0.4.0-beta.3

### Feat
- **arch:** add arch_imp.h interface header file which can be used in different architectures
- **arch:** change accelerate calculation interface to adapt newest code
- **config:** remove creation.create_t_set term (system temperature) in config file, use rescale instead.

### BREAKING CHANGE

creation.create_t_set term (system temperature) was removed in config file.


<a name="v0.4.0-beta.2"></a>
## [v0.4.0-beta.2] - 2019-10-05
### Build
- **dependency:** add yaml-cpp 0.6.2 package.

### Docs
- **changelog:** add changelog for version 0.4.0-beta.2.
- **readme:** update building and running document in readme.

### Feat
- **config:** add stages parsing in yaml config file.
- **config:** add config parsing for rescaling feature.
- **config:** remove config.toml and use config.yaml as config file.
- **config:** add a yaml formated config file: config.yaml.
- **dependency:** update kiwi lib to version 0.2.0.
- **random:** make LEGACY random a random number engines.
- **random:** add more random number generators (lcg, mt, stc, xoshiro, legacy).
- **random:** log random seed when using std::random_device to generate seed in auto mode.
- **rescale:** rescale system temperature by config in config file.

### Fix
- **collision:** fix bug of setting atom velocity in collision (wrong sqrt).
- **config:** fix config parsing of alloy ratio of Ni.

### Refactor
- **atom-list:** rename AtomList::size() to AtomList::cap().
- **config:** use yaml as config file type and add yaml config parsing.
- **configuration:** move common code in kinetic energy calculation to function mvv and function reduceEnergy.
- **rescale:** move resacle implementation in class WorldBuilder to namespace configuration.
- **simulation:** use stages to set variable-step-length and collision velocity.

### Test
- **collision:** add test for atom velocity setting in collision.

### BREAKING CHANGE

config file has been changed to yaml style.


<a name="v0.4.0-beta"></a>
## [v0.4.0-beta] - 2019-09-24
### Build
- **cmake:** add cmake config for building frontend unit tests.
- **libcomm:** update libcomm to v0.3.0.

### Ci
- **gitlab-ci:** update cache key in .gitlab-ci.yml file to match libcomm version change.

### Docs
- **changelog:** add changelog for version 0.4.0-beta.

### Feat
- **config:** add config term for outputting thermodynamic information(like temperature).
- **config:** raise error during config parsing if some required terms are not set.
- **config:** add config parsing for variable step length.
- **simulation:** set variable step length during simulation and log physical time.
- **thermodynamics:** output thermodynamic information in release mode which is controlled by config file.

### Fix
- **$compile:** fix compiling error "error: ‘memcpy’ is not a member of ‘std’" in gcc-7.2.
- **config:** fix wrong reference of config value passed to MDSimulation class.
- **logs:** add a space between colon and execution time in logs.
- **mem-leak:** fix warning of "deleting object that has virtual functions but non-virtual destructor

### Refactor
- **config:** move atoms dumping config to [output.dump] section in config file.
- **config:** use vector not pointer as type of members vsl_break_points and vsl_lengths in ConfigValues.
- **config:** move ConfigValues to frontend.
- **packer:** better log for unpacking received data in pack/LatPackerFirst class if size does n

### Test
- **config:** add tests for config parsing of variable step length.
- **neighbour:** add neighbor offset test to test indexing of neighbor atoms.

### BREAKING CHANGE

config file changed: new config term [output.thermo] must be specified in config file.

config file changed: atoms dumping config are moved to [output.dump] section.

config file changed: rename timesteps_length to def_timesteps_length, add
variable_step_length array.


<a name="v0.3.4"></a>
## [v0.3.4] - 2019-08-05
### Docs
- **changelog:** add changelog for version 0.3.4.

### Fix
- **output:** fix bug of NaN of atom position by setting block size and buffer size the same.

### Refactor
- **config:** move output configuration to struct Output, and set it as a member of class Config
- **output:** use a enum type OutputMode as output mode type.
- **output:** split debug dump and copy dump mode.
- **output:** move code of md outputting to frontend/io directory.
- **output:** move MDSimulation::output to output interface.
- **output:** move outputting simulation results to frontend.


<a name="v0.3.3"></a>
## [v0.3.3] - 2019-07-15
### Build
- **cmake:** modify sunway architecture configuration in cmake script.
- **sunway:** add cmake script to compile and link sunway slave code.
- **toolchain:** add cmake toolchain file for sunway taihulight supercomputer.

### Ci
- **gitlab-ci:** add gitlab-ci config file to build code.

### Docs
- **changelog:** add changelog for version 0.3.3.

### Feat
- **libcomm:** upgrade libcomm version to 0.2.0.
- **output:** add output of global inter atoms count and global lattice atoms count(in development mode).
- **simulation:** add onForceSolved callback for simulation, and move forceChecking() to frontend.
- **simulation:** add beforeStep and postStep callback in simulation.
- **version:** move version config to version.cmake file.
- **ws:** add algorithm implementation of voronoy diagram to ws utils.

### Refactor
- **eam:** inlining func sendForce in class atom into func atom::computeEam.
- **eam:** make atoms eam computing functions(latRho, interRho, latDf, latForce and interForce) private.
- **force:** move system force output(development mode only) to class MDSimulation in frontend.
- **ws:** add a common macro VORONOY to calculate coordinate of nearest lattice for a atom.

### Test
- **ws:** correct tests due to changes of implementation of ws algorithm.


<a name="v0.3.2"></a>
## [v0.3.2] - 2019-07-15
### Build
- **libcomm:** update libcomm version to v0.1.0 to fix ghost size bug.
- **pot:** update dependency lib pot version from v0.1.0-alpha to v0.1.0-beta.

### Docs
- **changelog:** add changelog for version 0.3.2.

### Feat
- **configuration:** add feature of system kinetic energy and temperature calculation for simulation.
- **logs:** add istty check for colorful logs.

### Fix
- **$compile:** fix compiling error of calling comm::neiSendReceive in openmpi env.
- **eam:** add missed rho contribution of neighbour atoms in the same bucket.
- **force:** add the missing force contributing of neighbour atoms in the same bucket(hash table) for inter atoms.
- **neighbour:** search more neighbour lattices(cut_lattice + 1) when making neighbour indexes.
- **neighbour:** fix bug of "Atom's force grows to be very large suddenly."
- **ws:** fix bug of converting atom postion to lattice index.

### Refactor
- **atoms-creator:** inlining function WorldBuilder::dofCompute().
- **configuration:** move systemForce function in class AtomSet to namespace configuration.
- **eam:** remove unused func params Domain* in eam calculating.
- **motion:** remove addtional function calling and redundant intermediate variables.

### Test
- **configuration:** add test for system temperature calculating of configuration::temperature.

### Pull Requests
- Merge branch 'develop' into 'master'
- Merge branch 'fix-missed-inter-force' into 'develop'


<a name="v0.3.1"></a>
## [v0.3.1] - 2019-04-23
### Build
- **cmake:** configure MD_DEV_MODE macro in  config.cmake file.
- **dependency:** update kiwi lib to version 0.2.0-beta.

### Docs
- **changelog:** add changelog for version 0.3.1

### Feat
- **cli:** add 'version' to command line options.
- **force:** add system force output in development mode.
- **force:** add implementation for systemForce computing in dev mode.
- **inter:** add more assert details when atom position is not as expected in InterParticlePacker::s
- **inter:** add member function addInterAtom/removeInterAtom and addGhostAtom to class InterAtomList.
- **lattice-atom:** add checking function for out-of-simulation-box atoms in lattice lists.
- **libcomm:** update libcomm version to v0.1.0-beta
- **libcomm:** add communication and domain dependency lib: libcomm.
- **neighbour:** add new neighbour lattice index and index iterator as well as its tests.

### Fix
- **$compile:** add MD_DEV_MODE macro to simulation.cpp/.h file to fix compilation issue.
- **$compile:** add missing header files to fix compilation issue.
- **eam:** use full neighbour index to calculate force and eam rho between inter atoms, instead of half neighbour index.
- **eam:** use reverted communication(z->y->x, not x->y->z) for electron density(rho) and force exchanging.
- **force:** fix the force reset bug: not all lattice atoms' force are reset.
- **ghost:** fix ghost atoms communication of inter atoms.
- **inter:** set df,rho,f and v values to 0 when receiving a new inter atom (or inter ghost atom).
- **inter:** remove mismatched coordinate checking method for received ghost atoms.
- **inter:** clear ghost inter atoms before each simulation step.
- **ws:** fix incorrect calculating of atom lattice coordinate in y and z dimension.

### Perf
- **inter:** add HashTable for inter(and ghost) atoms to reduce the time complexity of inter atoms searching.

### Refactor
- **eam:** use new neighbour atoms index to compute eam rho(electron density).
- **eam:** spilit atom::computeEam code into latRho,interRho,latDf,latForce,interForce.
- **force:** rename variable atom_neighbour_up to lattice_neighbour in inter force calculation.
- **force:** refactor func print_force to use lambda expressions.
- **force:** use new neighbour atoms index to compute force.
- **frontend:** move non md lib files to frontend directory.
- **lattice:** move members in AtomList to class BccLattice.
- **libcomm:** refactor AtomList::exchangeAtomFirst and exchangeAtom to use libcomm Packer.
- **libcomm:** remove getatomx/y/z funcs and use comm::fwCommLocalRegion to get lattice region for forwarding communication.
- **libcomm:** refactor InterAtomList::borderInter to use libcomm Packer.
- **libcomm:** refactor InterAtomList::exchangeInter to use libcomm Packer.
- **libcomm:** remove Crystal MD domain module and use domain module in libcomm.
- **packer:** mv pack and unpack funcs under namespace pack to corresponding Packer class.
- **packer:** remove temp variables in packer onSend and onReceive funcs in pre inlining commit.
- **packer:** simplify code of force unpack in ForcePacker::onReceive, LatPacker::onReceive and RhoPacker::onReceive.

### Test
- **ws:** add missing tests for ws_utils.


<a name="v0.3.0"></a>
## [v0.3.0] - 2019-03-11
### Build
- **$compile:** add debug in cmake config.
- **dependency:** move dependency to dependency.cmake.
- **googletest:** update googletest version to 1.8.1
- **pkg:** correct pkg dependency version of pot lib.
- **pkg:** update to pkg version 0.2.0-beta.
- **pkg:** use pkg tool version 0.2.0-alpha.

### Docs
- **README:** update authors scope in README.md
- **changelog:** add changelog for version 0.3.0
- **config:** add commnets for config parser.

### Feat
- **alloy:** create different atom types in atom_types.h
- **atom:** change basic type of sendlist and recvlist in domain from int to _type_atom_id.
- **atom:** add lambda feature to traverse all atoms in sub-box.
- **atom:** Merge branch 'feature-atom-element' into 'dev'
- **atom:** create a new struct AtomElement to store atom information like force,location,velocity, etc.
- **atom:** move atoms 1-d array to 3-d array.
- **atom-dump:** add atoms_dump_interval in config.
- **atom-dump:** add feature of dumping inter atoms in copy mode.
- **atom-dump:** add feature of dumping by frame.
- **atom-dump:** addd feature of dumping atoms before collision and dumping inter atoms to the same binary file.
- **atom-type:** add INVALID atom type to tag the atom that is out of lattice.
- **atom_creator:** add zreo momentum for multiple types of atoms.
- **changelog:** add changelog support
- **collision:** change the unit of PKA to eV.
- **config:** add config parsing for Fe-Cu-Ni alloy.
- **config:** add time_step_length configuration in config file.
- **datatype:** add pre_define.h file
- **denpendency:** use newest kiwi framework (bundle put/get without MPI_Pack & MPI_Unpack).
- **dependency:** update to kiwi framework.
- **domain:** add localBuild for domain.
- **domain:** add local boundary in domain.
- **eam:** add full eam adaptation for multi-atom-types.
- **eam:** use class EamPhiList to replace InterpolationObject array, and add better api.
- **eam:** using OneWayEamList to replace (InterpolationObject *rho) and (InterpolationObject *f).
- **inter:** add out-of-box check for inter atoms.
- **logs:** add kiwi logs support.
- **logs:** add logs to file support.
- **merge:** Merge branch 'feature-config-resolver' to branch 'master'
- **multi_atom_type:** using atom_type::atom_type as atom type.
- **potential:** add libpot and change code to use libpot api.
- **simulation-output test tools:** merge output files into one file using mpi-IO
- **test:** add test framework catch2 for testing.
- **tools:** conv tool can convert allmost binary atom file to text file.
- **ws:** add feature of obtaining lattice coordinate and index of nearest atom.
- **ws:** add isInBox member function for ws utils.

### Fix
- fix merge conflict of mergeing branch feature-multi-types-atom into dev.
- **$compile:** fix compile error.
- **atom:** fixed runtime segmentfault and incorrect results.
- **atom-dump:** fixed atom dump segment-fault by updating to newest kiwi.
- **atom-dump:** fixed bug: header of local storage is set many times.
- **atom_creator:** fixed bug [#3](https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/issues/3): duplicate atom id.
- **atoms:** fixed incrrect atoms count of system in atom::decide implementation.
- **collision:** fixed possible outofRange error in setv.
- **collision:** fixed segment error while packing particledata in exchangeInter.
- **compile:** fix compile error: pre_config.h not found.
- **config:** fixed simulation error, timesteplength not sync to other processors.
- **config:** fixed bug of hasError==true in clang compiler.
- **domain:** double the lattice size and lattice coord in x direction.
- **eam:** make (pointer to pointer) type to (pointer) type in eamBcast in eam.cpp.
- **eam:** fix index out-of-range bug in atom::computeEam when iterating the inter atoms.
- **inter:** we still add the unexcepted atom to inter atoms list InterAtomList::unpackExInterRecv to fix incorrct atom count issue.
- **inter:** fix incorrect iterator of InterAtomList::inter_list in atom.cpp
- **inter:** fix index-out-of-range bug of inter atoms in function computeEam due to incorrect periodic boundary.
- **inter:** fixed incrrect atoms count of system in InterAtomList::exchangeInter.
- **newton_motion:** remove fixed mass(only Fe mass is used) in newton_motion.cpp/.h
- **simulation:** filter the INVALID atoms when performing simulation.
- **test-domain:** fix compilation errors in tests of domain.
- **tools:** add filter for empty data in binary atom dump file.
- **tools:** fixed bugs of conv tool for processing binary file (result not right).

### Refactor
- rename src/config to src/types.
- remove some unnecessary code.
- **$config-resolver:** add mpi data pack util to resolve and synchronize configure informatio
- **atom:** move atom related files to "atom" directory.
- **atom:** move atom set(including atomlist and interatomlist) to class AtomSet.
- **atom-dump:** move atom dumping in atom.cpp to atom_dump.cpp
- **atom-list:** move atom-list relative mothods to class AtomList.
- **atom-set:** remode pointer of Domain in class AtomSet constructor.
- **atoms:** move member sendForce, sendrho, sendDfEmbed in class Domain to class atom.
- **atoms:** rename the index variable of neighbour atoms.
- **atoms_creater:** add a world builder (atoms builder) based on the old create_atom.
- **atoms_creator:** refactor computeScalar method in world_builder.cpp
- **decomposition:** rename variale in DomainDecomposition unified prefix.
- **domain:** add double-x lattice size and coord member for class Domain.
- **domain:** rename class DomainDecomposition to Domain.
- **domain:** move variable (nlocal nghost lolocal loghost) * (x,y,z) to DomainDecomposition.
- **domain:** move variable nlocalx, nlocaly, nlocalz in atom to domain.
- **domain:** extract domain boundary to class Region.
- **domain:** move code of sub-box bound in class atom to class DomainDecomposition.
- **domain:** refator doamin implementation to use lattice priority strategy (rather than measured length priority).
- **domain:** add Domain::Builder and move domain values to a class with all const values.
- **domaindecomposition:** add comments for domaindecomposition
- **eam:** extract eam parser to eam_parser.cpp/.h file.
- **eam:** move potential code to libpot.
- **inter:** use new added out-of-box checking method to exchange inter atoms with neighbour pro
- **inter:** move inter atoms to file inter_atom_list.h
- **inter:** move implementation of out-of-box inter atoms communication in InterAtomList to another
- **inter:** refactor class InterAtomList (atom/inter_atom.h/.cpp)
- **inter:** move more inter-atom-list relative methods into class InterAtomList.
- **inter:** change type of variable intersendlist and interrecvlist in class InterAtomList to type std::vector<std::vector<AtomElement *> >.
- **newton_motion:** move newton motion in atom.cpp/.h to newton_motion.cpp/.h
- **pack:** move pack and unzpck into "pack" directory.
- **pack:** refactor function unpack_recvfirst in pack.cpp
- **potential:** move potential file parsing in simulation.cpp to eam.cpp
- **potential:** separate potential file parsing for master processor and non-master processors.
- **simulation:** better code(comments, MACRO, etc.).
- **simulation:** refactor code of offset used for periodic boundary in exchangeAtom in domain.cpp
- **test:** switch from catch test framework to googletest.
- **ws:** move ws related functions to namespace ws.

### Test
- replace MPI_Init() with kiwi::mpUtils.
- **atom_creator:** add test for rescalar in WorldBuilder.
- **collision:** add MPI tests for particledata.
- **domain:** add test for domain decomposition
- **domain:** modified test code in domain_test.cpp
- **eam:** add test for phi(pair potentials).
- **ws:** add test of ws::isOutBox for periodic boundary.

### Pull Requests
- Merge branch 'dev' into 'master'

### BREAKING CHANGE

config Term collision_v have been removed.


<a name="v0.2.0"></a>
## [v0.2.0] - 2018-03-08
### Build
- **CMakeLists.txt pre_config.h:** add cmake install feature, add sunway support enable option.

### Docs
- **README.md .input.swo:** add contributing section to README.md and remove .input.swo file.


<a name="v0.1.3-sunway"></a>
## [v0.1.3-sunway] - 2017-11-17

<a name="v0.1.2-sunway"></a>
## [v0.1.2-sunway] - 2017-11-17

<a name="v0.1.1-sunway"></a>
## [v0.1.1-sunway] - 2017-11-17

<a name="v0.1.0"></a>
## v0.1.0 - 2017-11-17

[Unreleased]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.4.0...HEAD
[v0.4.0]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.4.0-beta.4...v0.4.0
[v0.4.0-beta.4]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.4.0-beta.3...v0.4.0-beta.4
[v0.4.0-beta.3]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.4.0-beta.2...v0.4.0-beta.3
[v0.4.0-beta.2]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.4.0-beta...v0.4.0-beta.2
[v0.4.0-beta]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.3.4...v0.4.0-beta
[v0.3.4]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.3.3...v0.3.4
[v0.3.3]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.3.2...v0.3.3
[v0.3.2]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.3.1...v0.3.2
[v0.3.1]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.3.0...v0.3.1
[v0.3.0]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.2.0...v0.3.0
[v0.2.0]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.3-sunway...v0.2.0
[v0.1.3-sunway]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.2-sunway...v0.1.3-sunway
[v0.1.2-sunway]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.1-sunway...v0.1.2-sunway
[v0.1.1-sunway]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.0...v0.1.1-sunway

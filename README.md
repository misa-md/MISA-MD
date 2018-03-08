# Crystal-MD

A molecular dynamics (MD) simulation program.  
Author:[Baihe](mailto:baihe_ustb@163.com)
Update:[Chugenshen](mailto:genshenchu@gmail.com)

### Build
#### Build from CMake (recommend)  
dependency:
1. cmake CMake 3.6 or higher is required;
2. c++ 11 feature required(check your gcc version);
3. mpi.

```bash
$ mkdir build
$ cd build
$ cmake ../
$ make

# insatll to path {install_prefix}, the default install path is /usr/local. 
# Use -DCMAKE_INSTALL_PREFIX={install_prefix} to change to another location.
$ make install
$ cp {install_prefix}/bin/CrystalMD ../example/  # copy executable file to example directory.
```

#### Build for sunway system
```bash
$ mkdir build
$ cd build
$ cmake ../ -DCMAKE_INSTALL_PREFIX={install_prefix} # default insatll path is /usr/local.
$ make
$ make install # insatll to position {install_prefix}.
$ cp {install_prefix}/bin/CrystalMD ../example/  # copy executable file to example directory.
```

#### Build from Makefile
```bash
$ mkdir build
$ make
$ cp ./CrystalMD ../example/  # copy executable file to example directory.
```

### Run
NOtice: befor running, you should have [FeCuNi.eam.alloy](https://www.ctcms.nist.gov/potentials/Download/Fe-Cu-Ni-GB/FeCuNi.eam.alloy) file,then specific the file path in config.toml file.  
For example, you can get the file by running following command to download the file:
```bash
$ wget https://www.ctcms.nist.gov/potentials/Download/Fe-Cu-Ni-GB/FeCuNi.eam.alloy -O exmaple/FeCuNi.eam.alloy
```

run simulation:
```bash
$ cd  example
$ ./CrystalMD  --help # run for showing help.
$ mpiexec -n 64 ./CrystalMD -c config.toml  # run MD simulation
```

### Contributing
It is meaningful to make commit messages formated, so we use use [AngularJS's commit message convention](https://github.com/angular/angular/blob/master/CONTRIBUTING.md#-commit-message-guidelines) also known as conventional-changelog.  
You can also use [commitizen tool](https://github.com/commitizen/cz-cli) to generate AngularJS style commit messages.
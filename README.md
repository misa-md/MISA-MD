# MISA-MD

A molecular dynamics (MD) simulation program.
Developers:[Baihe](mailto:baihe_ustb@163.com) and [Chugenshen](mailto:genshenchu@gmail.com)

![https://hpcer.pages.hpcer.dev/CrystalMD/MDoc/](https://img.shields.io/badge/document-online-ff69b4)

## Build
### Build from CMake (recommend)  
dependency:
1. CMake 3.6 or higher is required;
2. c++ 11 feature required(check your gcc version);
3. mpi.

we use [pkg](https://github.com/genshen/pkg) tool to manage c/c++ lib dependencies,
you must install pkg on your system.
To install dependencies in project root directory, run:
```
$ pkg fetch
$ pkg install
```

To build or install MISA-MD, run:
```bash
$ cmake -DCMAKE_BUILD_TYPE=Release -H. -Bbuild/
$ cmake --build build

# to install to path {install_prefix}, use -DCMAKE_INSTALL_PREFIX={install_prefix} to change to another location.
$ cmake --build build --target install
```

## Run
Notice: before running, you should have [FeCuNi.eam.alloy](https://www.ctcms.nist.gov/potentials/Download/Fe-Cu-Ni-GB/FeCuNi.eam.alloy) file,then specific the file path in config.toml file.  
see [here](https://www.ctcms.nist.gov/potentials/Fe-Cu-Ni.html) for more details.
For example, you can get the file by running following command to download the file:
```bash
$ wget https://www.ctcms.nist.gov/potentials/Download/Fe-Cu-Ni-GB/FeCuNi.eam.alloy -O exmaple/FeCuNi.eam.alloy
```

run simulation:
```bash
$ cd  example
$ ../build/bin/misamd  --help # run for showing help.
$ mpiexec -n 64 ../build/bin/misamd -c config.yaml  # run MD simulation
```

## Document
For details on how to build, configure and run, check out our [document](https://hpcer.pages.hpcer.dev/CrystalMD/MDoc/).

## Contributing
It is meaningful to make commit messages formatted, so we use [AngularJS's commit message convention](https://github.com/angular/angular/blob/master/CONTRIBUTING.md#-commit-message-guidelines) also known as conventional-changelog.  
You can also use [commitizen tool](https://github.com/commitizen/cz-cli) to generate AngularJS style commit messages.
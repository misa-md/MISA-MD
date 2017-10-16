# Crystal-MD

A molecular dynamics (MD) simulation program.  
Author:[Baihe](mailto:baihe_ustb@163.com)


### how to build&run

build from CMake (recommend)
```sh
mkdir build
cd build
cmake ../
make
mpiexec -n 100 ./Crystal-MD
```

build from Makefile
```sh
 make
 mpiexec -n 100 ./Crystal-MD
```

#!/bin/sh
set -e

# motivation: when using mpiswgcc to compile slave code via `-mslave` flag on sw5,
# error will raise: `Do not use -mhost and -mslave together`.
# Thus, in this script, it will replace mpiswgcc compiler with sw5gcc when there is a `-mslave` flag.
SWcc=sw5gcc
MPISWcc=mpiswgcc # mpi C compiler

declare -a finalopts
finalopts=()
is_slave=false

for o in "$@"; do
    if [ "$o" = "-mslave" ] ; then
        is_slave=true
    fi
    #add all other options to the list
    finalopts+=("$o")
done

if [ "$is_slave" = true ] ; then
  exec $SWcc "${finalopts[@]}"
else
    exec "$MPISWcc" "${finalopts[@]}"
fi

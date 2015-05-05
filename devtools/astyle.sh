#!/bin/sh

: ${BINDIR:="$( cd "$( dirname "$0" )" && pwd )"}

dir="${BINDIR}/../src"

exec astyle -r --options="${BINDIR}/astylerc" \
    --exclude="${dir}/utils" --exclude="${dir}/contrib" \
    "${dir}/*.h"  "${dir}/*.xmh" "${dir}/*.cc" "${dir}/*.cpp" "${dir}/*.c"  

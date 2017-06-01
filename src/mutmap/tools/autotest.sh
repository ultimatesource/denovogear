#!/bin/bash

# Client for automatically rebuilding when source files change. Start
# test_server.R first.

# This should be run from build directory like this:
# src/mutmap/autotest.sh

find ../src/mutmap/ -not -name "*.swp" -and -not -type d | entr sh -c 'ninja; curl localhost:6011'

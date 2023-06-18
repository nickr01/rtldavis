#!/bin/sh

./build.sh
$GOPATH/bin/rtldavis -gain 0 -maxmissed 21 -tf AU -noafc -dump


#!/bin/sh
set -e

git submodule init
git submodule update

go install -v .

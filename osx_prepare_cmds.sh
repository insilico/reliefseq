#!/bin/bash
#
# osx_prepare_cmds.sh - Bill White - 3/30/15
#
# source this command before 'configure" to set OSX-related paths'

export BOOST_ROOT="/usr/local"
export CXXFLAGS="-I/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/gcc/x86_64-apple-darwin13.4.0/4.9.2/include $CXXFLAGS"
export LDFLAGS="-L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9 $LDFLAGS"

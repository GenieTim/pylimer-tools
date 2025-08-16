#!/usr/bin/env bash

awk -v version="$1" '/## v/ {printit = substr($2, 1) == version}; /## main/ {printit = version == "main"}; printit;' "$2"

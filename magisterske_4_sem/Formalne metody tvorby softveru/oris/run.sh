#!/bin/bash

INSTALL_DIR=$( cd "$( dirname "$0" )" && pwd )
(cd "$INSTALL_DIR" && java -jar oris-runner-*.jar)
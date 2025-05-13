#!/bin/bash

# Modifications by S. Shepard
# Changes made for IRMA are on the IRMA@v1.0 branch

# For Linux.
git clone https://github.com/sammysheep/Complete-Striped-Smith-Waterman-Library.git \
    && cd Complete-Striped-Smith-Waterman-Library/src/ \
    && git checkout IRMA@v1.0 \
    && make \
    && ls -l ssw_v0.1.5M_x86_64

# If on Mac rename with Darwin64 suffix.
if [ "$(uname -s)" == "Darwin" ]; then
    mv ssw_v0.1.5M_x86_64 ssw_v0.1.5M_darwin64
fi

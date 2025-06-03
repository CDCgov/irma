#!/bin/bash

default_irma_core=v0.5.1

version=${PINNED_CORE:-$default_irma_core}
archive=irma-core-integrated-${version}.zip
url=https://github.com/CDCgov/irma-core/releases/download/${version}/$archive

cd IRMA_RES/scripts \
    && wget "$url" \
    && unzip "$archive" \
    && rm "$archive"

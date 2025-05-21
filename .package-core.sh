#!/bin/bash
version=$PINNED_CORE
archive=irma-core-integrated-${version}.zip
url=https://github.com/CDCgov/irma-core/releases/download/${version}/$archive

cd IRMA_RES/scripts \
    && wget "$url" \
    && unzip "$archive" \
    && rm "$archive"

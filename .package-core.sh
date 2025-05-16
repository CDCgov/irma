#!/bin/bash
version=${PINNED_CORE:-v0.5.0}
archive=irma-core-${version}-integrated.zip
url=https://github.com/CDCgov/irma-core/releases/download/${version}/$archive

cd IRMA_RES/scripts \
    && wget "$url" \
    && unzip "$archive" \
    && rm "$archive"

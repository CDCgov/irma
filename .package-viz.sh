#!/bin/bash

default_irma_viz=v0.2.0

version=${PINNED_VIZ:-$default_irma_viz}
archive=irma-viz-integrated-${version}.zip
url=https://github.com/CDCgov/irma-viz/releases/download/${version}/$archive

download_archive() {
    if command -v curl > /dev/null 2>&1; then
        curl -L --proto '=https' --tlsv1.2 -sSf --output "$archive" "$url"
    elif command -v wget > /dev/null 2>&1; then
        wget --https-only --secure-protocol=TLSv1_2 -q -O "$archive" "$url"
    else
        echo "Error: curl or wget is required to download: '$url'" >&2
        return 1
    fi
}

cd IRMA_RES/scripts \
    && download_archive \
    && unzip -o "$archive" \
    && rm "$archive"

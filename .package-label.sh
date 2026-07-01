#!/bin/bash

default_label_version=v0.7.2

version=${PINNED_LABEL:-$default_label_version}
archive=label-${version}-universal.zip
url=https://github.com/CDCgov/label/releases/download/${version}/$archive

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

[ -d "LABEL_RES" ] && rm -rf LABEL_RES LABEL

download_archive \
    && unzip "$archive" \
    && rm -rf flu-amd/LABEL_RES/training_data/{H3,H5,B_,H1,H7,H9}* flu-amd/LABEL_RES/scripts/creation \
    && mv flu-amd/{LABEL,LABEL_RES} . \
    && rm -rf flu-amd/ "$archive"

#!/bin/bash

version=${PINNED_LABEL:-v0.7.0}
archive=label-${version}-universal.zip
url=https://github.com/CDCgov/label/releases/download/${version}/$archive

[ -d "LABEL_RES" ] && rm -rf LABEL_RES LABEL

wget "$url" \
    && unzip "$archive" \
    && rm -rf flu-amd/LABEL_RES/training_data/{H3,H5,B_,H1,H7,H9}* flu-amd/LABEL_RES/scripts/creation \
    && mv flu-amd/{LABEL,LABEL_RES} . \
    && rm -rf flu-amd/ "$archive"

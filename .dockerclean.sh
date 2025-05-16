#!/bin/bash

TP=IRMA_RES/third_party
kill_arch=aarch64
if [[ "$(uname -m)" != "x86_64" ]]; then
    kill_arch=x86_64
fi

rm $TP/*Darwin $TP/*$kill_arch || exit 1

if [ -d "LABEL_RES/third_party" ]; then
    TP=LABEL_RES/third_party
    rm $TP/*Darwin $TP/*$kill_arch || exit 1
fi

rm IRMA_RES/scripts/irma-core_*{Darwin,$kill_arch}

#!/bin/bash
bash3.2 -eo pipefail \
    IRMA FLU /samples/SRR38262444_1.fastq /samples/SRR38262444_2.fastq "$(basename "$0" .sh)" > "$0.log" 2>&1 \
    && echo -e "Test success:\t$0" \
    || {
        cat "$0.log"
        echo -e "Test failure:\t$0"
        exit 1
    }

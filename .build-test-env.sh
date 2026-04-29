#!/bin/bash
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    --file Dockerfile.test-env \
    --tag ghcr.io/cdcgov/irma/test-env:latest \
    --push \
    .

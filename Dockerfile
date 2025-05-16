FROM debian:bookworm-slim AS builder


RUN apt-get update --allow-releaseinfo-change --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y wget git unzip zip ca-certificates \
    && update-ca-certificates


ENV PINNED_LABEL=v0.7.0
ENV PINNED_CORE=v0.5.0

ARG core_branch

WORKDIR /app
COPY . .
RUN ./.package-label.sh \
    && ./.package-core.sh \
    && ./.dockerclean.sh \
    && rm ./.*.sh

FROM debian:bookworm-slim AS base

ARG APT_MIRROR_NAME=
RUN if [ -n "$APT_MIRROR_NAME" ]; then sed -i.bak -E '/security/! s^https?://.+?/(debian|ubuntu)^http://'"$APT_MIRROR_NAME"'/\1^' /etc/apt/sources.list && grep '^deb' /etc/apt/sources.list; fi
RUN apt-get update --allow-releaseinfo-change --fix-missing \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y procps zip perl \
    && apt clean autoclean \
    && apt autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

WORKDIR /app
COPY  --from=builder /app /app

ENV PATH="/app:${PATH}"
WORKDIR /data
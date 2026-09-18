#!/bin/bash
# Useful for installing needed, local resources after `git pull`.
# Not meant to be used for packaging for deployment.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

os=$(uname -s)
arch=$(uname -m)

case "$os" in
    Darwin)
        # Apple binaries are universal (fat) and already cover x86_64 + arm64
        keep_suffix="Darwin"
        remove_suffixes=(Linux_x86_64 Linux_aarch64)
        ;;
    Linux)
        case "$arch" in
            x86_64) other_arch=aarch64 ;;
            aarch64 | arm64) other_arch=x86_64 ;;
            *)
                echo "INSTALL_LOCAL ERROR: unsupported architecture '$arch'" >&2
                exit 1
                ;;
        esac
        keep_suffix="Linux_$arch"
        remove_suffixes=(Darwin "Linux_$other_arch")
        ;;
    *)
        echo "INSTALL_LOCAL ERROR: unsupported OS '$os'" >&2
        exit 1
        ;;
esac

echo "==> Installing LABEL"
./.package-label.sh

echo "==> Installing IRMA-core"
./.package-core.sh

echo "==> Installing IRMA-viz"
./.package-viz.sh

echo "==> Removing binaries for other platforms ($os/$arch keeps '$keep_suffix')"
# Only remove unusable binaries that are not committed.
shopt -s nullglob
[ -d LABEL_RES/third_party ] && for suffix in "${remove_suffixes[@]}"; do
    rm -f LABEL_RES/third_party/*"$suffix"
done

for suffix in "${remove_suffixes[@]}"; do
    rm -f IRMA_RES/scripts/irma-core_*"$suffix" IRMA_RES/scripts/irma-viz_*"$suffix"
done

echo "==> Install complete."

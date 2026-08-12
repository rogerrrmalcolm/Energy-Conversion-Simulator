#!/bin/sh
set -eu

if ! grep -qw avx2 /proc/cpuinfo; then
    echo "error: this image requires an x86-64 host with AVX2 support" >&2
    exit 78
fi

exec /opt/black-hole/bin/black-hole-sim "$@"

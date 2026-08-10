# syntax=docker/dockerfile:1

FROM ubuntu:24.04 AS builder

ARG BH_ENABLE_AVX2=ON

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        cmake \
        ninja-build \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
COPY . .

RUN cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DBH_BUILD_TESTS=OFF \
        -DBH_ENABLE_AVX2=${BH_ENABLE_AVX2} \
    && cmake --build build --parallel \
    && cmake --install build --prefix /opt/black-hole

FROM ubuntu:24.04 AS runtime

COPY --from=builder /opt/black-hole/bin /opt/black-hole/bin
COPY --from=builder /opt/black-hole/share /opt/black-hole/share

LABEL org.opencontainers.image.title="Black Hole Energy Simulation" \
      org.opencontainers.image.description="C++20 Kerr and Penrose-process simulation CLI"

USER 65532:65532
WORKDIR /opt/black-hole

ENTRYPOINT ["/opt/black-hole/bin/black_hole_demo"]
CMD ["--help"]

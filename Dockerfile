# Stage 1: Build Stage
FROM ubuntu:22.04 AS builder

# Install dependencies and remove install files to reduce image size
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Set working directory in the container for build
WORKDIR /ffea_src

# Copy the project files into the container
COPY . .

# Create a build directory and configure/build/install ffea using CMake
RUN mkdir -p build && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/ffea/ \
    && cmake --build . --parallel 4 \
    && cmake --install .

# Stage 2: Final Runtime Image
FROM ubuntu:22.04

# Install runtime dependencies only
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Copy the installed files from the builder stage
COPY --from=builder /ffea /ffea

# Add the two binary directories to path
# So that ffea can be used anywhere, and default python is venv
ENV PATH="/ffea/bin:/ffea/venv/bin:${PATH}"
ENV VIRTUAL_ENV="/ffea/venv"

# Set working directory for user
WORKDIR /workspace

# Define default command (optional, can be overridden)
#CMD ["ffea"]
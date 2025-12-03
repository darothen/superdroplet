#!/bin/bash
#
# Quick build script for sd_cpp
#

set -e

echo "Building sd_cpp..."

# Create build directory if it doesn't exist
if [ ! -d "build" ]; then
    mkdir build
fi

cd build

# Configure with CMake
echo "Configuring..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building..."
cmake --build . -j4

echo ""
echo "Build complete! Run with: ./build/sd_cpp"
echo ""

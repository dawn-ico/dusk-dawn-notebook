#!/usr/bin/env sh

# deployment of AtlasUtils isn't clear yet, so we just install a trimmed
# down version here:

mkdir build

g++ \
  -c AtlasExtractSubmesh.cpp \
  -o build/AtlasExtractSubmesh.o \
  -isystem /usr/local/atlas/include \
  -isystem /usr/local/eckit/include \
  -std=c++17

g++ \
  -c AtlasCartesianWrapper.cpp \
  -o build/AtlasCartesianWrapper.o \
  -isystem /usr/local/atlas/include \
  -isystem /usr/local/eckit/include \
  -std=c++17

g++ \
  -c GenerateRectAtlasMesh.cpp \
  -o build/GenerateRectAtlasMesh.o \
  -isystem /usr/local/atlas/include \
  -isystem /usr/local/eckit/include \
  -std=c++17

ar \
  rcs build/libatlasUtilsLib.a \
  build/AtlasExtractSubmesh.o \
  build/AtlasCartesianWrapper.o \
  build/GenerateRectAtlasMesh.o

mkdir -p /usr/local/AtlasUtils/lib
cp build/libatlasUtilsLib.a /usr/local/AtlasUtils/lib/
mkdir -p /usr/local/AtlasUtils/include/atlas_utils/utils
cp *.h /usr/local/AtlasUtils/include/atlas_utils/utils/

# Building with CMake

## build with Clang

- Dependence:
    - cmake version 3.31.3 or newer
    - Ninja
    - Clang
    - OneapiTBB (for parallel)

### Debug
```sh
cmake --preset="clang-debug" --no-warn-unused-cli
cmake --build build/clang/Debug
```

### Release
```sh
cmake --preset "clang-release" --no-warn-unused-cli
cmake --build build/clang/Release
```


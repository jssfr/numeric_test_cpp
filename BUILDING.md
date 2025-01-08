# Building with CMake

## build with Clang

- Depends:
    - cmake version 3.31.3 or newer
    - Ninja
    - Clang

### Debug
```sh
cmake --preset="clang-debug" --no-warn-unused-cli
cmake --build build/clang/Debug
```

Here are the steps for building in release mode with a multi-configuration
generator, like the Visual Studio ones:

### Release
```sh
cmake --preset "clang-release" --no-warn-unused-cli
cmake --build build/clang/Release
```


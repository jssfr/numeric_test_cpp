# Building with CMake

## build with Clang

- Dependence:
    - cmake version 4.0.2 or newer
    - Ninja
    - Clang, Gcc or msvc
    - Gcc 15.1 or newer
    - Clang 20.1 or newer
    - OneapiTBB (for parallel) (optinal)

### Debug
```sh
cmake --preset clang-debug --no-warn-unused-cli
cmake --build --preset bclang-debug
```

```sh
cmake --preset gcc-debug --no-warn-unused-cli
cmake --build --preset bgcc-debug
```

```sh
cmake --preset msvc --no-warn-unused-cli
cmake --build --preset bmsvc --config Debug
```


### Release
```sh
cmake --preset clang-release --no-warn-unused-cli
cmake --build --preset bclang-release
```

```sh
cmake --preset gcc-release --no-warn-unused-cli
cmake --build --preset bgcc-release
```

```sh
cmake --preset msvc --no-warn-unused-cli
cmake --build --preset bmsvc --config Release
```

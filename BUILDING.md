# Building with CMake

## build with Clang

- Dependence:
    - cmake version 3.31.3 or newer
    - Ninja
    - Clang or msvc
    - OneapiTBB (for parallel)

### Debug
```sh
cmake --preset clang-debug --no-warn-unused-cli
cmake --build --preset bclang-debug
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
cmake --preset msvc --no-warn-unused-cli
cmake --build --preset bmsvc --config Release
```

# Ray Tracing in One Weekend - C

With my renewed interest in C, I decided to re-do [Ray Tracing in One Weekend](https://raytracing.github.io/) in C. I might also use OpenGL/OpenCL in the future.

Other Ray Tracing in One Weekend repos:

- https://github.com/gau-nernst/ray-tracing-numba
- https://github.com/gau-nernst/ray-tracing-rust

On Linux and macOS

```bash
make main
make main ENABLE_OPENMP=1  # for OpenMP support
./main 0
```

On MacOS, Apple Clang does not support OpenMP. Install LLVM from Homebrew to build with OpenMP.

```bash
brew install llvm
make main CC=$(brew --prefix llvm)/bin/clang ENABLE_OPENMP=1
```

On Windows

Download [`stb_image.h`](https://github.com/nothings/stb/blob/master/stb_image.h) and place it under `include/external`.

```bash
cl /Iinclude /Zc:preprocessor /std:c11 /O2 /fp:fast /Femain.exe src/*.c
cl /Iinclude /Zc:preprocessor /std:c11 /O2 /fp:fast /openmp /Femain.exe src/*.c  # for OpenMP support
./main.exe 0
```

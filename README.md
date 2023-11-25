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

## Learnings

- Implement `Vec3`:
  - Container: having `Vec3` be a union (C11) of a struct and an array makes it convenient to access individual fields explicitly or through a loop. Accessing inactive fields (i.e. fill values via struct members but access values via array) is allowed in C (no undefined behavior).
  - Operations: since C does not support operator overloading, we have to implement each operation as a seperate function, including different combinations of operands (Vec3 + Vec3 or Vec3 + float), and explicitly call individual functions when we do calculations. The code can look very verbose, but we can improve it in two ways:
    1. Thanks to C11's `_Generic()`, we can define a macro to select an appropriate function based on operands type.
    2. We can use variadic macro magic (with `__VA_ARGS__`) to make + operation accept more than two operands. I got the solution from [this Stackoverflow answer](https://stackoverflow.com/a/11763277). Sadly, since C does not support recursive macro, we need to pre-define the maximum number of operands a `vec3_add()` macro can take. In my case, supporting up to 4 operands is enough.
- Pseudo-random generator (PRNG): I went to this repo [lemire/testingRNG](https://github.com/lemire/testingRNG) to look for a simple and fast PRNG. Since I'm using `float` instead of `double`, I only need a 32-bit PRNG. Moreover, since MSVC does not have a native 128-bit unsigned integer type (`__uint128_t` in GCC and Clang), I would need to write platform-conditioned code (MSVC has some compiler intrinsics to deal with 128-bit integers like `_mul128()`)
- Object-oriented programming (OOP):
  - Since the raytracing in one weekend series adopts an OOP style, I need to implement interface in C (in particular, for `Hittable`, `Material`, and `Texture`). From my research on the Internet, there are 3 main ways to do this:
    1. Add an enum tag as the first field of a struct, and class-specific data as subsequent fields. An interface's method will check this tag to call the appropriate function (using a `switch` statement).
    2. Store function pointers as the first few fields of a struct. There must be an initialization function to assign these function pointers correctly.
    3. Store address of a vtable in a struct. The vtable contains addresses of all methods of a class. This is similar to approach 2, but now each struct only needs to store 1 address. This is only beneficial if the interface has a lot of methods.
  - In any approaches, we need to cast pointer type from interface pointer type to a particular class pointer type. This is unsafe since there is no guarantee that data in that address is from that valid class. If we need some elements of OOP, it's probably better to use a programming language that supports it natively i.e. C++ or Rust.
  - I'm not sure if it's possible to write a useful (and simple to use) raytracer without OOP. We can have a list of primitive types (e.g. `Sphere`) and loop over it. However, with OOP, many features can be effortlessly implemented by referencing another object implementing an interface (e.g. `Translate`, or `Checker`).

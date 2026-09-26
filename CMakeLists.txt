cmake_minimum_required(VERSION 4.3)
project(Cpp_Prototype)

set(CMAKE_CXX_STANDARD 20)

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)

# optimization flags
add_compile_options(
        -O3               # Maximum speed optimization
        -march=native     # Target your specific Intel Skylake/Coffee Lake CPU architecture
        -mavx2            # 256-bit wide SIMD vector registers (8 floats per cycle)
        -mfma             # Fused multiply-add in single clock cycle
        -fblocks          # Required for Apple Grand Central Dispatch (dispatch_apply)
        -ffast-math       # Enables aggressive floating-point vectorization
)

add_compile_options(-Rpass=loop-vectorize)
add_executable(Cpp_Prototype main.cpp
        Matrix.h
        AdamsBashforth.h
        AdamsBashforth.cpp
        CrankNicolson.h
        CrankNicolson.cpp
        Chorin.h
        Chorin.cpp
        Renderer.cpp
        Renderer.h
        Utils.h
        Utils.cpp
        StreamFunction.h
        StreamFunction.cpp
        Obstacle.h
        Obstacle.cpp
        ImmersedBoundary.h
        ImmersedBoundary.cpp
)

target_compile_options(Cpp_Prototype PRIVATE -fsanitize=address -g)
target_link_options(Cpp_Prototype PRIVATE -fsanitize=address)

find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
add_compile_definitions(ACCELERATE_NEW_LAPACK)

add_library(glad glad.c)
target_include_directories(glad PUBLIC external/glad/include)

target_link_libraries(${PROJECT_NAME} PRIVATE
        "-framework Accelerate"
)

target_link_libraries(${PROJECT_NAME}
        PRIVATE
        glad
        glfw
        OpenGL::GL
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
)

cmake_minimum_required(VERSION 3.20)
project(vma1)

add_compile_options(
        "$<$<CONFIG:RELEASE>:-O3>"
        "$<$<CONFIG:DEBUG>:-O0>"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_INCLUDE_CURRENT_DIR true)

# ThreadSanitizer
set(TSAN_FLAGS
        "-fsanitize=thread"
        CACHE STRING "Flags used by the C++ compiler during ThreadSanitizer builds."
        FORCE)

# AddressSanitize
set(ASAN_FLAGS
        "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer"
        CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
        FORCE)

# LeakSanitizer
set(LSAN_FLAGS
        "-fsanitize=leak -fno-omit-frame-pointer"
        CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
        FORCE)

# MemorySanitizer
set(MSAN_FLAGS
        "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer"
        CACHE STRING "Flags used by the C++ compiler during MemorySanitizer builds."
        FORCE)

# UndefinedBehaviour
set(UBSAN_FLAGS
        "-fsanitize=undefined"
        CACHE STRING "Flags used by the C++ compiler during UndefinedBehaviourSanitizer builds."
        FORCE)

# Build Types
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${ASAN_FLAGS}")

add_executable(vma1 main.cpp
        TimeMeasurer/time_measurer.cpp
        TimeMeasurer/time_measurer.h
        Algebra/algebra.cpp
        Algebra/algebra.h Matrix/abstract_matrix.h Matrix/diagonal_box_matrix.h Matrix/mutable_matrix.h Matrix/transposed_matrix.h Algebra/gmres.h)

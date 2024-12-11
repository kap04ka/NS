cmake_minimum_required(VERSION 3.10)

# Название проекта
project(NS)

# Устанавливаем стандарт C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Указание директорий для заголовков
include_directories(${CMAKE_SOURCE_DIR}/include)

# Рекурсивно собираем все файлы из src/
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Создаем исполняемый файл, добавляя только SOURCES (HEADERS можно не добавлять)
add_executable(NS ${SOURCES})

# Проверяем наличие OpenMP
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    # Добавляем флаги компиляции и связывания для OpenMP
    target_link_libraries(NS PUBLIC OpenMP::OpenMP_CXX)
endif()

# Опционально: добавляем флаги компиляции, если нужны
target_compile_options(NS PRIVATE -Wall -Wextra -O2)

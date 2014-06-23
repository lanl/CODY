add_executable(square
               main.cpp
               host.cpp
               device.cpp
               app.hpp
               ${CMAKE_SOURCE_DIR}/src/common/app-base.hpp)


add_test(square square -p -d -v)

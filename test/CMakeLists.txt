#test/CMakeLists.txt

INCLUDE_DIRECTORIES( include ${CMAKE_SOURCE_DIR}/include )

ADD_EXECUTABLE(plate plate.cpp)
ADD_EXECUTABLE(slab slab.cpp)
ADD_EXECUTABLE(slab2 slab2.cpp)
ADD_EXECUTABLE(slab3 slab3.cpp)


TARGET_LINK_LIBRARIES(plate mySrc)
TARGET_LINK_LIBRARIES(slab mySrc)
TARGET_LINK_LIBRARIES(slab2 mySrc)
TARGET_LINK_LIBRARIES(slab3 mySrc)
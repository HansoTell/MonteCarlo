# alle cpp Dateien automatisch sammeln
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

add_library(MonteCarloLib ${LIB_SOURCES})

target_include_directories(MonteCarloLib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

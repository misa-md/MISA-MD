
add_library(${MD_PLUGINS_LIB_NAME} STATIC
        plugin_api.h
)

target_include_directories(
        ${MD_PLUGINS_LIB_NAME}
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
        ${MD_PLUGINS_LIB_NAME}
        PUBLIC
        ${MD_LIB_NAME}
)

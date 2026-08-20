if(CMAKE_CONFIGURATION_TYPES)
    set(BIN_DEST_3RD "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}")
else()
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(BIN_DEST_3RD "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}")
    else()
        set(BIN_DEST_3RD "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}")
    endif()
endif()

set(COPY_COMMANDS)
foreach(SOURCE_DIR IN LISTS BIN_3RD_LIST)
    list(APPEND COPY_COMMANDS
        COMMAND ${CMAKE_COMMAND} -E make_directory "${BIN_DEST_3RD}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SOURCE_DIR}" "${BIN_DEST_3RD}"
    )
endforeach()

add_custom_target(_ready_bin
    ${COPY_COMMANDS}
    COMMENT "Copying OCCT runtime DLLs next to executables"
    VERBATIM
)

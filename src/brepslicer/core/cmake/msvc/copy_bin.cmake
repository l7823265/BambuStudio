# Multi-config generators: stage OCCT DLLs next to every configuration output.
if(CMAKE_CONFIGURATION_TYPES)
    set(_OCC_DESTS
        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}"
        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}"
    )
else()
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(_OCC_DESTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}")
    else()
        set(_OCC_DESTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}")
    endif()
endif()

set(COPY_COMMANDS)
foreach(SOURCE_DIR IN LISTS BIN_3RD_LIST)
    foreach(DEST IN LISTS _OCC_DESTS)
        list(APPEND COPY_COMMANDS
            COMMAND ${CMAKE_COMMAND} -E make_directory "${DEST}"
            COMMAND ${CMAKE_COMMAND} -E copy_directory "${SOURCE_DIR}" "${DEST}"
        )
    endforeach()
endforeach()

add_custom_target(_ready_bin
    ${COPY_COMMANDS}
    COMMENT "Copying OCCT runtime DLLs next to executables"
    VERBATIM
)

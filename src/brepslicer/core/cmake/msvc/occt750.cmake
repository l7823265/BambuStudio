# Open CASCADE 7.5.0 — same install used by D:/OCC (VS2019 / vc14).
# This machine currently provides Debug binaries only (libd / bind).

set(OCCT750_DIR "D:/lib/opencascade-7.5.0/build_vs2019/win64" CACHE PATH "OCCT 7.5 win64 root")

if(NOT EXISTS "${OCCT750_DIR}/inc/Standard.hxx")
    message(FATAL_ERROR "OCCT headers not found at ${OCCT750_DIR}/inc")
endif()

set(OCCT750_INCLUDE_DIR "${OCCT750_DIR}/inc")
set(OCCT750_LIBRARY_DIR "${OCCT750_DIR}/vc14/libd")
set(OCCT750_BINARY_DIR  "${OCCT750_DIR}/vc14/bind")

if(NOT EXISTS "${OCCT750_LIBRARY_DIR}/TKernel.lib")
    message(FATAL_ERROR "OCCT debug libraries not found at ${OCCT750_LIBRARY_DIR}")
endif()

file(GLOB OCCT750_LIBS "${OCCT750_LIBRARY_DIR}/*.lib")

add_library(OCCT750 INTERFACE)
target_include_directories(OCCT750 INTERFACE ${OCCT750_INCLUDE_DIR})
target_link_directories(OCCT750 INTERFACE ${OCCT750_LIBRARY_DIR})
target_link_libraries(OCCT750 INTERFACE ${OCCT750_LIBS})

list(APPEND BIN_3RD_LIST ${OCCT750_BINARY_DIR})

list(LENGTH OCCT750_LIBS _occt_lib_count)
message(STATUS "OCCT 7.5: ${OCCT750_DIR}")
message(STATUS "OCCT include: ${OCCT750_INCLUDE_DIR}")
message(STATUS "OCCT libs:    ${OCCT750_LIBRARY_DIR} (${_occt_lib_count} .lib)")

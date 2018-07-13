
file(GLOB sources *.f90)
list(REMOVE_ITEM sources "main.f90")

add_library(libtreecode ${sources})
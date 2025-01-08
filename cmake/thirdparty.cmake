# This commands downloads AND configures fmt. It sets up some variables.
CPMAddPackage(
    NAME fmt
    GIT_TAG 11.1.1
    VERSION 11.1.1
    GITHUB_REPOSITORY fmtlib/fmt
    SOURCE_DIR ${THIRDY_DIR}/fmt
)
#CPMAddPackage("gh:fmtlib/fmt#11.1.1") # other way of use


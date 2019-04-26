@REM Creates Visual Studio solution using CMake
@REM
@REM Syntax: create_devenv arg1
@REM
@REM arg1  : defines the CMake generator that will be used to create the
@REM         solution files. Please check the supported generators with
@REM         CMake -h.
@REM
@REM Notes:
@REM
@REM The command %~1 removes the quotes from arg1, as it can contain white
@REM space characters. This is needed because otherwise there would be
@REM double quotes when commands are expanded by the shell.
@REM
@REM The build directory is derived from the given generator by replacing
@REM the white space characters with underscore. The syntax for Windows
@REM command shell is %X: =_% where X is an environment variable.
@REM

@REM Identify positional arguments
@set CMAKE_GENERATOR=%~1
@set SOURCE_DIR=%~dp0
@set BUILD_DIR=%~dp0_build\%CMAKE_GENERATOR: =_%

@if exist "%BUILD_DIR%" (
    @echo Cleaning build directory: %BUILD_DIR%
    rmdir /S /Q "%BUILD_DIR%"
)

@mkdir "%BUILD_DIR%"
@pushd "%BUILD_DIR%"
cmake -G "%CMAKE_GENERATOR%" "%SOURCE_DIR%"

@if %ERRORLEVEL% NEQ 0 (
    @echo.
    @echo Something went wrong. Please check the logs.
    @echo.
    exit /B 1
)

start Teisko.sln

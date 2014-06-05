@echo OFF
ECHO.
ECHO This batch file is used to invoke the QC script by drag'n'dropping of
ECHO a MaxQuant txt folder (or any file within) onto this .bat file.
ECHO Within the same txt folder, a QC report in PDF format will be created.
ECHO.
ECHO First time installation of this script [for admins only]:
ECHO The script expects a certain structure to be present.
ECHO Right next to this .bat file should by a subfolder named '_internal'
ECHO which holds:
ECHO      'R-3.1.0'             containing a complete R installation 
ECHO                            (including all packages required to run the 
ECHO                             PTXQC package)
ECHO      'compute_QC_report.R' The R script that is called by this .bat
ECHO.
REM configure this directory according to above instructions
REM 
REM  --- config done --- 
REM

REM A usage message with multiline (i.e. keep the empty line!)
setlocal EnableDelayedExpansion
set USAGE= ^

Usage: ^

%0 ^<txt-folder^> [config.yaml] ^

Please try again.
REM end of USAGE

REM check number of arguments (must be 1 [txt folder] or 2 [txt folder + yaml])
set argC=0
for %%x in (%*) do Set /A argC+=1
if %argC%==0 (
  ECHO Not enough arguments^^! (No arguments provided^)
  ECHO !USAGE!
  goto end
)
if !argC! gtr 2 (
  ECHO Too many arguments! (At most two are allowed^)
  ECHO !USAGE!
  goto end
)



set I=%~dp0\_internal
REM use drive+path of this .bat file 
ECHO.
ECHO Found QC directory at %I%
ECHO.

REM Check upon first argument (the txt folder - or any file in it)
if not exist %1 ( 
  goto not_found
) 
set is_dir=0
REM PushD+PopD is the only reliable way I found to check if %1 is a directory or a file
REM Everything else proposed (e.g., %1\. or %1\nul does not work on Win7)
PUSHD %1 && POPD || goto is_file
set is_dir=1
:is_file
if %is_dir%==1 ( 
  echo -- %1 is a directory 
  set txt=%1
) else ( 
  echo -- %1 is a file
  set txt=%~dp1
) 

ECHO Txt folder is at '%txt%'

if %argC%==2 (
  set yaml_file=%2
  ECHO.
  REM check file extension. We could also use %~x2 to get the file extension
  if !yaml_file:~-5! NEQ .yaml (
    ECHO Error: Second argument '!yaml_file!' does not seem to be a '.yaml' file.
    ECHO Please rename or use proper file.
    goto end
  )
  if exist !yaml_file! (
    ECHO Using YAML configuration file !yaml_file!
    ECHO.
  ) else (
    ECHO Error: Yaml file '!yaml_file!' given as second argument does not exist.
    ECHO Every QC script run will create a default .yaml file in the respective txt folder, which you can COPY
    ECHO and use as input to any subsequent run.
    ECHO.
    goto end
  ))

REM use R_LIBS, not R_LIBS_USER, since the latter will APPEND to the search path, i.e. if the package is installed locally, it will take precedence and mess
REM up versioning
set R_LIBS=%I%\R-3.1.0\library
REM echo R_LIBS=%I%\R-3.1.0\library

%I%\R-3.1.0\bin\x64\rscript --vanilla %I%\compute_QC_report.R !txt! !yaml_file!
REM report error, if any
if ERRORLEVEL 1 (
  ECHO Failed folder %txt%
  ECHO Failed folder %txt% >> \error.log
)

goto end

:not_found
echo Could not find the given (network) folder '%1'^^! Please contact your admin/bioinformatician^^!

:end
pause;


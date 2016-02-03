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
  ECHO Too many arguments^^! (At most two are allowed^)
  ECHO Args:
  for %%x in (%*) do echo %%x
  ECHO !USAGE!
  goto end
)

REM usually, path without quotes are needed, e.g. for set PATH=!I!;%PATH%  (spaces are no problem)
set I=%~dp0\_internal
REM sometimes quoted strings are required (i.e. when calling rscript)
set Iqt="%~dp0\_internal"

REM currently, we do not allow Drag'n'drop PTXQC to be installed on a path with spaces, since
REM this breaks pandoc.exe (and the temp directories it creates) for some reason
for /f "tokens=2" %%a IN ("%I%") do goto NoSpaces

REM use drive+path of this .bat file 
ECHO.
ECHO Found QC directory at '%I%'
ECHO.

REM Check upon first argument (the txt folder - or any file in it)
REM use bare %1 not "%1", it works fine, even with spaces
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
  echo -- '%1' is a directory 
  set txt=%1
) else ( 
  echo -- '%1' is a file
  set txt=%~dp1\
) 

ECHO Txt folder is at '%txt%'

if %argC%==2 (
  REM we use %~dpnx2 instead of %2 to get rid of potential extra surrounding quotes
  set yaml_file="%~dpnx2"
  ECHO.
  REM check file extension. We could also use %~x2 to get the file extension
  if !yaml_file:~-5! NEQ .yaml (
    ECHO Error: Second argument '!yaml_file!' does not seem to be a '.yaml' file.
    ECHO Please rename or use proper file.
    goto end
  )
  if exist !yaml_file! (
    ECHO Using YAML configuration file '!yaml_file!'
    ECHO.
  ) else (
    ECHO Error: Yaml file '!yaml_file!' given as second argument does not exist.
    ECHO Every QC script run will create a default .yaml file in the respective txt folder, which you can COPY
    ECHO and use as input to any subsequent run.
    ECHO.
    goto end
  )
)

REM use R_LIBS, not R_LIBS_USER, since the latter will APPEND to the search path, i.e. if the package is installed locally, it will take precedence and mess
REM up versioning
set R_LIBS=!I!\R-3.1.0\library
REM echo R_LIBS=%I%\R-3.1.0\library

set PATH=!I!\R-3.1.0\bin\pandoc;!PATH!

where pandoc
if ERRORLEVEL 1 (
  ECHO.
  ECHO 'Pandoc.exe not found in '!I!\R-3.1.0\bin\pandoc\pandoc.exe'.
  ECHO.
  ECHO Please install Pandoc from 'https://github.com/jgm/pandoc/releases'
  ECHO and copy it to the folder mentioned above.
  ECHO.
  goto end
)

REM works WITH and WITHOUT spaces in !txt!
if "%yaml_file%" NEQ "" (
  ECHO Using YAML '!yaml_file!' file and calling R now ...
)
  ECHO.
  ECHO Calling '%I%\R-3.1.0\bin\x64\rscript.exe --vanilla %I%\compute_QC_report.R !txt! !yaml_file!'
  REM using Iqt to guard against path with spaces. Note that using a manually quoted !I!, i.e. "!I!", does not work
  CMD /C "!Iqt!\R-3.1.0\bin\x64\rscript --vanilla !Iqt!\compute_QC_report.R !txt! !yaml_file!"
)

set myscriptEL=%errorlevel%
REM report error, if any
if %myscriptEL% GTR 0 (
  ECHO Failed folder %txt%. Does it exist and is writeable^^!?
  ECHO Failed folder %txt%. Does it exist and is writeable^^!? >> \error.log
  ECHO.
  ECHO The error should be indicated above.
  ECHO Common errors are ^(yours might be different^)
  ECHO   - rscript.exe is not present in the folder mentioned above
  ECHO   - your R installation is not 64bit and runs out of memory
  ECHO   - the file '%I%\compute_QC_report.R' is missing
  ECHO.
  ECHO If the error persists, please make a screenshot of this window and the folders mentioned above and open an issue ticket
  ECHO on https://github.com/cbielow/PTXQC. 
  ECHO You can also copy the text from this window by 1^) right-clicking here, 2^) choose 'Select all', 3^) press the 'Enter' key, 4^) paste it somewhere ^(issue ticket, email^).
  ECHO This ensures that all text is available to debug the problem.
  ECHO Thanks^^!
)

goto end

:noSpaces
ECHO Drag'n'Drop PTXQC is installed in a path containing spaces ("!I!")
ECHO This is not allowed (since it breaks HTML report generation). Please move it to a path without spaces.
ECHO.
goto end

:not_found
echo Could not find the given (network) folder '%1'^^! Please contact your admin/bioinformatician^^!
goto end

:end
pause






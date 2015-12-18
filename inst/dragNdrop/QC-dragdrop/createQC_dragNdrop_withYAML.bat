@echo OFF
REM Here we hardcode the YAML file which is forwarded to 'createQC_dragNdrop.bat'
REM The '%~dp0' just prefixes the YAML file with the path to this .bat file
REM i.e. the YAML file should be right next to it.
REM You can also put the .yaml file somewhere else and hardcode the path
REM e.g. set yaml_file=c:\temp\my.yaml
REM Just make sure the file ending remains '.yaml'
set yaml_file=%~dp0\config.yaml

ECHO.
ECHO This batch file allows invoking 'createQC_dragNdrop.bat'
ECHO with a YAML config file which resides in this directory.

REM check number of arguments (must be exactly 1 [txt folder])
set argC=0
for %%x in (%*) do Set /A argC+=1
if %argC% NEQ 1 (
  ECHO.
  ECHO Wrong number of arguments^!
  ECHO.
  ECHO Exactly one argument (containing the txt folder or any file within^) is expected^!
  ECHO.
  goto end
)


ECHO.
ECHO Using '%yaml_file%' as YAML config
ECHO.
ECHO.
REM Calling original bat file with two arguments
REM We use %~dp0 to get the path of this bat file, since the working directory depends on the calling environment
"%~dp0\createQC_dragNdrop.bat" %1 "%yaml_file%"


:end
pause;  
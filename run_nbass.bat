SET NBASS_PATH="%HOMEPATH%\nbass_distro"
SET NBASS_RUNS_PATH="%NBASS_PATH%\runs"
SET IDL_PATH="%NBASS_PATH%:%NBASS_PATH%/lib:%NBASS_PATH%/tools:%IDL_PATH%"
idl run_nbass.pro -quiet

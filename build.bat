REM Building script for the SenoCETSA project
REM ensure an appropriate conda environment is active before running

pyinstaller --name SenoCETSA^
  --add-data "README.md:."^
  --add-data "LICENSE:."^
  --add-data "cetsa_config.toml:."^
  --hide-console minimize-early^
  cetsa_interface.py
  
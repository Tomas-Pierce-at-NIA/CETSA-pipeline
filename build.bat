REM Building script for the SenoCETSA project
REM ensure an appropriate conda environment is active before running

pyinstaller --name SenoCETSA^
  --add-data "README.md:."^
  --add-data "LICENSE:."^
  --add-data "cetsa_config.toml:."^
  --add-data "icons:icons"^
  --hide-console minimize-early^
  --icon "icons/SenoCETSAicon_notext.ico"^
  cetsa_interface.py
  
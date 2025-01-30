REM Building script for the SenoCETSA project

conda activate senocetsa

pyinstaller --name SenoCETSA^
  --add-data "README.md:."^
  --add-data "LICENSE:."^
  --add-data "cetsa_config.toml:."^
  --hide-console minimize-early^
  cetsa_interface.py
  
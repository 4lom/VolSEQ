# VolSEQ

## Introduction

## How to use
Launch the python script using:

```
python3 VolSEQ.py
```

## Packaging

### Windows

To create a package pyinstaller is needed.

Launch the BuildExe.bat in power shell:

```
.\BuildExe.bat VolSEQ.py
```

### Linux

To create a package pyinstaller is needed.

Launch the BuildExe.bat in a shell:

```
.\BuildExe.sh VolSEQ.py
```

### MacOS

To create a package pyinstaller is needed.

Launch the BuildExe.bat in a shell:

```
.\BuildExe.sh VolSEQ.py
```

## Requirement for developement

python >=3.9 (not tested on older version)

Import needed:
 - Base64
 - PySimpleGUI
 - pathlib
 - io
 - PIL
 - pydicom
 - matplotlib

To download the needed import use pip:

1. Install PIL

```
pip install pillow
```

2. Install PIL

```
pip install PySimpleGUI
```

3. Install pydicom

```
pip install pydicom
```

4. Install matplotlib

```
python -m pip install -U matplotlib
```

5. Install pathlib

```
pip install pathlib
```
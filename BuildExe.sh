#!/bin/sh
pyinstaller -F --onefile --clean --windowed --collect-submodules=pydicom %1
# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = []
hiddenimports += collect_submodules('pydicom')


block_cipher = None


a = Analysis(
    ['VolSEQ.py'],
    pathex=[],
    binaries=[],
    datas=[ ('LICENSE', '.') 
            ],
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

app = BUNDLE( exe,
            name='VolSEQ.app',
            icon=None,
            bundle_identifier=None
)


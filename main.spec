# -*- mode: python -*-

block_cipher = None


a = Analysis(['bem.py'],
             pathex=[],
             binaries=[],
             datas=[("icon_bem.ico","."),
                    ("xfoil_executables/xfoil.exe","xfoil_executables"),
                    ("xfoil_executables/xfoil","xfoil_executables"),
                    ("karlsen.bem",".")],
             hiddenimports=['scipy._lib.messagestream','scipy.special.cython_special'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

Key = ['mkl']

def remove_from_list(input, keys):
    outlist = []
    for item in input:
        name, _, _ = item
        flag = 0
        for key_word in keys:
            if key_word in name:
                flag = 1
                break
        if flag != 1:
            outlist.append(item)
    return outlist

print("List of binaries:")
for _name,_path,_type in a.binaries:
    print(_name,_path,_type)

a.binaries = remove_from_list(a.binaries, Key)

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='main',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=False,
          runtime_tmpdir=None,
          console=False,
          icon="icon_bem.ico")

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=False,
               name='main')
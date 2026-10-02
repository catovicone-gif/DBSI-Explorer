# Builds dist\DBSI-Explorer.exe, a single-file Windows executable of the GUI
# (no Python installation needed on the machine that runs it).
#
# One-time setup, from the repo root:
#     python -m venv .venv
#     .venv\Scripts\python -m pip install -r requirements.txt pyinstaller
# Then:
#     powershell -ExecutionPolicy Bypass -File build_exe.ps1

$root = $PSScriptRoot
& "$root\.venv\Scripts\python.exe" -m PyInstaller `
    --noconfirm --clean --onefile --windowed `
    --name "DBSI-Explorer" `
    --paths "$root\DBSI_GUI" --paths "$root\bin" `
    --add-data "$root\DBSI_GUI\ArxivPaperR04InsolationExpansion.pdf;." `
    --exclude-module PyQt5 --exclude-module PyQt6 --exclude-module PySide2 --exclude-module PySide6 `
    --exclude-module IPython --exclude-module scipy --exclude-module pandas `
    --distpath "$root\dist" --workpath "$root\build" --specpath "$root\build" `
    "$root\DBSI_GUI\dbsi_gui_app.py"
